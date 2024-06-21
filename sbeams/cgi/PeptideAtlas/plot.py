#!/proteomics/sw/python/python3/bin/python3
#modified from https://analyticalsciencejournals.onlinelibrary.wiley.com/doi/10.1002/pmic.202100103
import cgi
import re, io,os
from io import BytesIO
import os.path
import sys
import subprocess
import json
from urllib.parse import parse_qs
os.environ['MPLCONFIGDIR'] = '/net/dblocal/www/html/devZS/sbeams/tmp/'

import pandas as pd
import numpy as np
from scipy.stats import linregress
import seaborn as sns
import statsmodels.stats.power as power
from tabulate import tabulate
import matplotlib.pyplot as plt
from sklearn.impute import KNNImputer
import json


sys.path.append("/net/dblocal/www/html/devZS/sbeams/lib/python/PeptideAtlas/extern")
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.io as pio
#import ext.easyMLR as emlr
#script_directory = os.path.dirname(os.path.abspath(__file__))


def equalize_medians(arr,log10=False):
    if log10:
        arr = np.log10(arr)
        arr = arr[np.isfinite(arr)]

    # Check if any column has NaN values
    if np.any(np.isnan(arr)):
        # Calculate the median while ignoring NaN values
        median_arr = np.nanmedian(arr, axis=0)
        # Find columns with NaN values
        nan_columns = np.isnan(median_arr)
        # Replace NaN median values with 0
        median_arr[nan_columns] = 0
        # get the max median
        median_global = np.max(median_arr)
    else:
        # Calculate the median of the entire array
        median_global = np.median(arr)

    # Find the difference between the median and the median of each column
    deltas = median_global - np.nanmedian(arr, axis=0)

    # Add the differences to each column
    arr_equalized = arr + deltas
    return arr_equalized

def plot_histograms(data=pd.DataFrame,log10=False):
    # transform data and set Series name
    if log10:
        data = np.log10(data)
        data = data[np.isfinite(data)]
    # Draw plot
    n=int((len(data.columns) +5)/5)*5

    fig, axes = plt.subplots(nrows=int(n/5), ncols=5,  figsize=(15, 15))
    axes = axes.ravel()
    idx=0
    for col in data:
        if idx == len(data.columns):
           break
        axes[idx].hist(data[col])
        axes[idx].set_title(data[col].name)
        idx +=1

    while idx < n: 
        fig.delaxes(axes[idx])
        idx +=1

    plt.tight_layout()
    return fig
def plot_boxplot(data=pd.DataFrame,log10=True):
    # transform data and set Series name
    #plt.figure(figsize=(8, 6))
    if log10:
        data = np.log10(data)
        data = data[np.isfinite(data)]
    #data = data.iloc[:, :5];
    #print(data.iloc[:, :5].to_string(index=False))
    try:
        fig=sns.boxplot(data=data)
        #fig.tick_params(axis='x', labelrotation=45)
        plt.xticks(rotation=45, ha='right',fontsize=10);
        return fig
    except Exception as e:
        return (f"An error occurred: {e}")

def plot_densityplot(data=pd.DataFrame,log10=True):
    # transform data and set Series name
    if log10:
        data = np.log10(data)
        data = data[np.isfinite(data)]
    #data = data.iloc[:, :5];
    #print(data.iloc[:, :5].to_string(index=False))


    rows_with_missing = data[data.isnull().any(axis=1)]
    rows_without_missing = data[~data.isnull().any(axis=1)]
    try:
        fig,axes=plt.subplots(1, 2, figsize=(14, 7), sharex=True, sharey=True)
        # Density plot for rows with missing values
        i=0
        for column in rows_with_missing.columns:
            density=sns.kdeplot(rows_with_missing[column].dropna(), ax=axes[0], label=column)

        x_peak = density.get_lines()[0].get_xdata()[np.argmax(density.get_lines()[0].get_ydata())]
        y_peak = np.max(density.get_lines()[0].get_ydata())
        # Add a vertical line at the peak
        axes[0].axvline(x=x_peak, color='black', linestyle='--', label=f'Peak: {x_peak:.2f}', linewidth=2)


        axes[0].set_title('Density Plot for Proteins with Missing Values')
        axes[0].set_xlabel('Intensity log(2)')
        axes[0].legend()
        # Density plot for rows without missing values
        for column in rows_without_missing.columns:
            density = sns.kdeplot(rows_without_missing[column], ax=axes[1], label=column)

        x_peak = density.get_lines()[0].get_xdata()[np.argmax(density.get_lines()[0].get_ydata())]
        y_peak = np.max(density.get_lines()[0].get_ydata())
        # Add a vertical line at the peak
        axes[1].axvline(x=x_peak, color='black', linestyle='--', label=f'Peak: {x_peak:.2f}', linewidth=2)

        axes[1].set_title('Density Plot for Proteins without Missing Values')
        axes[1].set_xlabel('Intensity (log2)')
        axes[1].legend()
        
        return fig
    except Exception as e:
        return (f"An error occurred: {e}")


def plot_sample_correlations_sns(df: pd.DataFrame,
                             data_columns: str = "I_.*",
                             correlation_function: callable = lambda x: np.corrcoef(x.T),
                             mode: str = "scatter", log10: bool = True):
    """
    Generates either a grid of paired scatter plots, or a heatmap of pairwise correlation values.
    
    Parameters
    ----------
    df: pandas DataFrame
    data_columns: regular expression string, default = "Intensity (.*)"
        Regular expression matching to columns containing data and a capture group for sample labels.
    correlation_function: callable, default = lambda x: np.corrcoef(x.T)
        Callable function to calculate correlation between samples.
    mode: str, default = "scatter"
        One of 'scatter' and 'heatmap' to swtich modes.
    log10: bool = True
        If True, data is log10 transformed prior to analysis
    
    Returns
    -------
    a plotly Figure
        either a subplotted, or a single heatmap depending on the mode
    """
    # pick and process data
    df_sub = df[[el for el in df.columns if re.match(data_columns, el)]].copy()
    if log10:
        df_sub = df_sub.apply(np.log10)

    df_sub = df_sub.replace([np.inf, -np.inf], np.nan)
    df_sub.columns = [re.findall(data_columns, el)[0] for el in df_sub.columns]

    if mode == "scatter":
        fig=sns.pairplot(df_sub, kind='reg', diag_kind='kde')
        # Iterate over the axes of the pair plot
        for i, ax in enumerate(fig.axes.flat):
            if i % (len(df_sub.columns) + 1) == 0:  # Skip the diagonal plots
                continue
            # Create a boolean mask for NaN values
            
            # Get the column names of the variables being plotted
            x_var = df_sub.columns[i // len(df_sub.columns)]
            y_var = df_sub.columns[i % len(df_sub.columns)]
            
            x_values=df_sub[x_var]
            y_values=df_sub[y_var]
            x_valid_indices = ~np.isnan(x_values)
            y_valid_indices = ~np.isnan(y_values)
            valid_indices = x_valid_indices & y_valid_indices
            x_valid = x_values[valid_indices]
            y_valid = y_values[valid_indices]
            # Calculate regression statistics
            slope, intercept, r_value, p_value, std_err = linregress(x_valid, y_valid)
            #print (f"slope={slope} intercept={intercept}" )
            # Annotate the plot with regression values
            annotation = f'Slope: {slope:.2f}\nIntercept: {intercept:.2f}\nR-value: {r_value:.2f}'
            ax.annotate(annotation, xy=(0.03, 0.85), xycoords='axes fraction', fontsize=8, color='red')
    elif mode=="heatmap":
        # Remove columns with all NaN values
        df_sub = df_sub[~np.isnan(df_sub).any(axis=1)]
        correlation_matrix = np.corrcoef(df_sub, rowvar=False)
        fig=sns.clustermap(correlation_matrix, cmap='coolwarm', annot=False, fmt=".2f", 
                                       xticklabels=list(df_sub.columns),yticklabels=list(df_sub.columns) )
        plt.xticks(fontsize=5, fontweight='bold', fontfamily='serif')
        plt.yticks(fontsize=5, fontweight='bold', fontfamily='serif')
        plt.setp(fig.ax_heatmap.xaxis.get_majorticklabels(), rotation=30)
        plt.setp(fig.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
        #plt.title('Correlation Heatmap', fontsize=14, fontweight='bold', fontfamily='serif')
        #fig.update_layout(coloraxis1=dict(showscale=False), font_size=5, width=800, height=800)
    else:
        raise ValueError
    return fig

def plot_sample_correlations(df: pd.DataFrame,
                             data_columns: str = "I_.*",
                             correlation_function: callable = lambda x: np.corrcoef(x.T),
                             mode: str = "scatter", log10: bool = True, binning: int = 10):
    """
    Generates either a grid of paired scatter plots, or a heatmap of pairwise correlation values.
    
    Parameters
    ----------
    df: pandas DataFrame
    data_columns: regular expression string, default = "Intensity (.*)"
        Regular expression matching to columns containing data and a capture group for sample labels.
    correlation_function: callable, default = lambda x: np.corrcoef(x.T)
        Callable function to calculate correlation between samples.
    mode: str, default = "scatter"
        One of 'scatter' and 'heatmap' to swtich modes.
    log10: bool = True
        If True, data is log10 transformed prior to analysis
    binning: int = 10
        Only relevant in scatter mode. Number of bins per full axis unit (e.g. 1e10 Intensity)
        used for datashader density rendering.
    
    Returns
    -------
    a plotly Figure
        either a subplotted, or a single heatmap depending on the mode
    """
    # pick and process data
    df_sub = df[[el for el in df.columns if re.match(data_columns, el)]].copy()
    if log10:
        df_sub = df_sub.apply(np.log10)

    df_sub = df_sub.replace([np.inf, -np.inf], np.nan)
    df_sub.columns = [re.findall(data_columns, el)[0] for el in df_sub.columns]

    if mode == "scatter":
        # setup subplots and axes
        fig = make_subplots(rows=len(df_sub.columns), cols=len(df_sub.columns), start_cell='bottom-left',
                            shared_yaxes=True, shared_xaxes=True, horizontal_spacing=0.03, vertical_spacing=0.03)
        i_range = (np.floor(np.nanmin(df_sub)), np.ceil(np.nanmax(df_sub))+1/binning)
        j_range = (np.floor(np.nanmin(df_sub)), np.ceil(np.nanmax(df_sub))+1/binning)
        i_width = int((i_range[1]-i_range[0]-1/binning)*binning+1)
        j_width = int((j_range[1]-j_range[0]-1/binning)*binning+1)
        print (f"i_range={i_range} j_range={j_range}")
        # fill plots
        for i,ni in enumerate(df_sub.columns):
            for j,nj in enumerate(df_sub.columns):
                # apply datashader
                dc = ds.Canvas(plot_width=i_width, plot_height=j_width, x_range=i_range, y_range=j_range)
                df_ij = df_sub[[ni,nj]].dropna() if i!=j else pd.DataFrame(df_sub[ni].dropna())
                da = dc.points(df_ij, x=ni, y=nj)
                zero_mask = da.values == 0
                da.values = da.values.astype(float)
                da.values[zero_mask] = np.nan
                
                # add trace
                fig.add_trace(
                    go.Heatmap(z=da,coloraxis="coloraxis1" if i!=j else "coloraxis2"),
                    row=j+1, col=i+1
                )
                
                # add annotations
                if j == 0:
                    fig.update_xaxes(title_text=ni, row=j+1, col=i+1, tickvals=list(range(0,i_width,binning)),
                                     ticktext=np.round(da[nj].values[0:i_width:binning]))
                if i == 0:
                    fig.update_yaxes(title_text=nj, row=j+1, col=i+1, tickvals=list(range(0,j_width,binning)),
                                     ticktext=np.round(da[ni].values[0:j_width:binning]))
                if i!=j:
                    fig.add_annotation(dict(text=str(np.round(np.min(correlation_function(df_sub[[ni,nj]].dropna())),4)),
                                            x=binning, y=j_width, showarrow=False), row=j+1, col=i+1)
        
        # layout figure
        fig.update_layout(template="simple_white", coloraxis2=dict(showscale=False, colorscale=["black", "black"]),
                          width=i*200+100, height=j*200+50, margin_t=0)
    elif mode=="heatmap":
        da = np.ones((len(df_sub.columns), len(df_sub.columns)))
        for i,ni in enumerate(df_sub.columns):
            for j,nj in enumerate(df_sub.columns):
                 #filter data and store correlation values
                df_ij = df_sub[[ni,nj]].dropna() if i!=j else pd.DataFrame(df_sub[ni].dropna())
                if i!=j:
                    da[i,j] = np.round(np.min(correlation_function(df_sub[[ni,nj]].dropna())),4)
        # create figure and label axes
        fig = go.Figure(data=go.Heatmap(z=da))
        fig.update_xaxes(tickvals=list(range(0,i+1,1)),
                          ticktext=list(df_sub.columns))
        fig.update_yaxes(tickvals=list(range(0,j+1,1)),
                          ticktext=list(df_sub.columns))
        fig.update_layout(template="simple_white", width=i*50+100, height=j*50+100)
    else:
        raise ValueError
    return fig
def plot_fc(df: pd.DataFrame,min_fc: float, fdr=0.01   ):

    # extract data
    #df =  df[df['Experiments'] == pat]
    df = df.rename (columns={'Biosequence Name': 'protein'})
    df.iloc[:, 3:] = df.iloc[:, 3:].astype(float)
    
    df.insert(df.shape[1], "square_cutoff",
               ["non_sig" if abs(fc)<min_fc or qval>fdr else "sig" for fc,qval in zip(df.fc, df.qval)])
    #sd = max([df_in[sc1].std(axis=1).mean(),
    #          df_in[sc2].std(axis=1).mean()])
    #exp_pow = power.tt_ind_solve_power(effect_size=min_fc/sd, alpha=fdr,
    #                                   nobs1=len(sc1))
    # Extract unique labels
    labels = df['Experiments'].unique()
    # Create subplots: the number of rows is the number of unique labels
    num_subplots = len(labels)
    num_cols = 3 if num_subplots > 3 else num_subplots
    num_rows = (num_subplots - 1) // num_cols + 1
    width = 800 if 400 * num_cols < 800 else 400*num_cols
    fig = make_subplots(rows=num_rows, cols=num_cols, subplot_titles=labels, shared_xaxes=True, shared_yaxes=True)

    # Generate scatter plots for each labels
    for i, label in enumerate(labels): 
        # Split the DataFrame based on labels
        df_sub = df[df['Experiments'] == label]
        row = i // num_cols + 1
        col = i % num_cols + 1 
        subplot=px.scatter(x=df_sub.fc,
                y=-np.log10(df_sub['pval'].astype(float)),
                color=df_sub.square_cutoff,
                template='simple_white',
                hover_name=df_sub.protein,
                #labels=dict(x="log2 fold change", y="-log10(p-value)",
                #            color="FDR {}, power {}".format(str(fdr), str(np.round(exp_pow, 2)))
                #           )
                ).update_layout(title={'text':f"<b>{label}</b>",
                                           'x':0.5,
                                           'xanchor': 'center',
                                           'y': 0.9,
                                           'yanchor': 'top'
                                           },
                                    title_font=dict(size=12),
                                    margin=dict(l=50, r=50, t=100, b=100))
        for trace in subplot.data:
            fig.add_trace(trace, row=row, col=col)

        fig.add_vline(min_fc, line_width=1, row=row, col=col)\
	    .add_vline(-min_fc, line_width=1, row=row, col=col)\
            .add_hline(-np.log10(df_sub.loc[df_sub.square_cutoff=="sig", "pval"].max()), line_width=1, row=row, col=col)

    fig.update_layout(showlegend=False)

    # Update layout to add titles
    fig.update_layout(
	xaxis_title='log2 fold change',
	yaxis_title='-log10(p-value)',
        height=600*num_rows,
        width=width,
        legend=dict(x=1.05, y=1),
        showlegend=True 
    )
    # Update x-axis title position to the middle of the subplot layout
    fig.update_xaxes(title_standoff=25)

    #svg_buffer = BytesIO()
    #fig.write_image(svg_buffer, format='svg')
    # Get the SVG content from the buffer
    #svg_content = svg_buffer.getvalue().decode('utf-8')

    #return svg_content
    return fig

def run_ttest(df: pd.DataFrame,
              c1_str: str, c2_str: str,
              cols_ann: list = ["Majority protein IDs", "Gene names", "Annotation"],
              s0: float or None = 0.05, fdr: float = 0.01, min_fc: float or None = None,
              n_perm: int = 2, plot_fdr_line: bool = True):
    """
    Run two-sided student's t-test with permutation based FDR control, either generating s0-based volcano lines or
    power-based quare cutoffs.
    
    Parameters
    ----------
    df: pandas DataFrame
        DataFrame containing data and annotations in pivot format
    c1, c2: regular expression string
        Regular expressions to extract columns containing replicate data for condition 1 and condition 2
    cols_ann: list, default = ["Majority protein IDs", "Gene names", "Annotation"]
        Additional annotation columns used for hover_annotation
    s0: float or None, default = 0.05
        If this is set and min_fc is None, non-linear volcano cutoff lines will be generated.
    fdr: float, default = 0.01
        False discovery rate for differentially expressed proteins.
    min_fc: float or None, default = None
        If the minimum fold change is set and s0 is None, square significance cutoffs will be used and the
        experimental power for the given fold change is calculated and reported.
    n_perm: int, default = 2
        Number of permutations to perform for FDR control. Recommended values are much higher.
    plot_fdr_line: bool, default = True
        Only effective for non-linear volcano lines. If true volcano lines are drawn. This is included because
        drawing these lines can takes some time.
    
    Returns
    -------
    a pandas DataFrame
        containing the selected data and annotation columns, and added columns for statistical output:
        fc, tval, pval, tval_s0, pval_s0, qval, FDR {fdr}%
    a plotly Scatter figure
        highlighting significant proteins and reporting FDR and s0 or power respectively.
    """
    c1 = re.compile(re.escape(c1_str))
    c2 = re.compile(re.escape(c2_str))

    # extract data
    df_in = df[[col for col in df.columns if col in cols_ann or re.search(c1, col) or re.search(c2, col)]].copy()
    ## remove rows with >=20% missing values
    proteins =df_in['protein']
    columns_not_protein = [col for col in df_in.columns if 'I_' in col]
    df_subset = df_in[columns_not_protein]
    threshold = 0.8  # 80% non-null values required
    df_cleaned_subset = df_subset.dropna(thresh=threshold*len(df_subset.columns))
    # Concatenate the first column back to the cleaned DataFrame
    df_cleaned = pd.concat([proteins[df_cleaned_subset.index], df_cleaned_subset], axis=1)
    df_in = df_cleaned

    ##replace nan value with column min
    for col_name  in df_in.columns:
        # Find the minimum value in the column ignoring NaN values
        if col_name == 'protein':
            continue
        
        min_value = np.nanmin(df_in[col_name])
        # Replace NaN values with the minimum value
        df_in[col_name].fillna(min_value, inplace=True)
        #df_in[col_name] = 10 ** df_in[col_name]

    df_in.reset_index(drop=True, inplace=True)
    # non-linear volcano lines
    sc1=[col for col in df_in.columns if re.search(c1, col)]
    sc2=[col for col in df_in.columns if re.search(c2, col)]

    
    # non-linear volcano lines
    if s0 and not min_fc:
        res, fig = emlr.perform_ttest_analysis(df_in, id_col=cols_ann[0],
                                               c1 = sc1,
                                               c2 = sc2, 
                                               plot_fdr_line=plot_fdr_line, s0=s0, fdr=fdr, n_perm=n_perm)
        fig.update_layout(legend_title_text="FDR {}, s0 {}".format(str(fdr), str(s0)))
    
    # square cutoffs (first performs ttest with s0 = 0, then adds power analysis)
    elif not s0 and min_fc:
        res, _ = emlr.perform_ttest_analysis(df_in, id_col=cols_ann[0],
                                             c1 = sc1, 
                                             c2 = sc2, 
                                             plot_fdr_line=False, s0=0, fdr=fdr, n_perm=n_perm)
        res.insert(res.shape[1], "square_cutoff",
                   ["non_sig" if abs(fc)<min_fc or qval>fdr else "sig" for fc,qval in zip(res.fc, res.qval)])
        sd = max([df_in[sc1].std(axis=1).mean(),
                  df_in[sc2].std(axis=1).mean()])
        exp_pow = power.tt_ind_solve_power(effect_size=min_fc/sd, alpha=fdr,
                                           nobs1=len(sc1))
        fig = px.scatter(x=res.fc,
                         y=-np.log10(res.pval),
                         color=res.square_cutoff,
                         template='simple_white', 
                         render_mode="svg",
                         hover_name=res.protein,
                         labels=dict(x="log2 fold change", y="-log10(p-value)",
                                     color="FDR {}, power {}".format(str(fdr), str(np.round(exp_pow, 2)))
                                    )
                        ).update_layout(width=600, height=700, 
                                        title={'text':f"<b>{c1_str} vs {c2_str}</b>",
                                               'x':0.5,
                                               'xanchor': 'center',
                                               'y': 0.9,  
                                               'yanchor': 'top'
                                               },
                                        title_font=dict(size=12),
                                        margin=dict(l=50, r=50, t=100, b=100))\
        .add_vline(min_fc, line_width=2).add_vline(-min_fc, line_width=2)\
        .add_hline(-np.log10(res.loc[res.square_cutoff=="sig", "pval"].max()), line_width=2)

    return res, fig
def sort_dataframe_result(result, column_name='fc', asc=False):
    try:
        sorted_result = result.sort_values(by=column_name, ascending=asc)
        return sorted_result
    except KeyError:
        print("Column", column_name, "not found in the DataFrame.\nresultset not sorted\n")
        return result
def print_svg(id):
    buffer = io.BytesIO()
    plt.savefig(buffer, format='svg', bbox_inches='tight')
    buffer.seek(0)
    # Read the SVG content
    svg_content = buffer.getvalue().decode()
    # Find the <svg> element and add the id attribute
    svg_content = svg_content.replace('<svg ', f'<svg id="{id}" ')
    print(svg_content)
def print_svg_plotly(fig):
    #svg_bytes = pio.to_image(fig, format='svg')
    #print(svg_bytes.decode())
    # Save the plot to a JSON string
    plot_json = fig.to_json()
    print (plot_json)

def read_data(file):
    # Execute Perl code using subprocess
    perl_code ='read_query_resultset';

    process = subprocess.Popen(['perl', perl_code, file],
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE, 
                                    text=True)
    # Check if the command was successful

    stdout, stderr = process.communicate()
    if process.returncode == 0:
        data = json.loads(stdout)
        # Extract column names from the first row
        columns = data[0]
        df = pd.DataFrame(data[1:],columns=columns)
        #print (df)
    else:
        print("Error executing Perl script.")
        #print("Error message:", result.stderr)
    return df

def knn_impute(df):
    # Convert data to numpy array for imputation, skipping the first row and first two columns
    data = df.iloc[1:,2:]

    # Replace '' with np.nan
    data = data.replace('', np.nan)

    # Perform KNN imputation
    imputer = KNNImputer(n_neighbors=3)
    imputed_data = imputer.fit_transform(data)


    # Combine imputed data with original first two columns
    # final_data = [data[0]]
    first_column_array = df['protein'].to_numpy()
    first_column_array = np.delete(first_column_array, 0, axis=0)
    second_column_array = df['Biosequence Gene Name'].to_numpy()
    second_column_array = np.delete(second_column_array, 0, axis=0)

    formatted_data = np.round(imputed_data, 4)
    imputed_data =  np.column_stack((first_column_array , second_column_array , formatted_data))
    return imputed_data


#### For command-line usage
def main():
    # HTML response

    ## test
    #set_name = 'protI_guest_20240531-130325-138';
    #query_params= {'mode': 'box_plot'}

    # Get query parameters from URL
    query_string = os.environ.get("QUERY_STRING", "")
    query_params = parse_qs(query_string)
    set_name = ''
    plot_type = ''

    # Parse form data
    form = cgi.FieldStorage()
    # Print input parameters
    #for key in form.keys():
    #    value = form.getvalue(key)
    #    print(f"{key}: {value}")
    if form is not None:
        if 'rs_set_name' in form.keys():
            set_name = form.getvalue('rs_set_name')
        if 'mode' in query_params:
            plot_type = query_params['mode']
    for idx, arg in enumerate(sys.argv[1:]):
        if idx==0:
            plot_type = arg
        elif idx==1:
            set_name = arg

    if set_name is None or plot_type is None:
        print (f"ERROR: failed to get the resultset file or plot type")
        return

    proteins = read_data(set_name)
    #print ( proteins.columns);

    data = proteins[[col for col in proteins.columns if re.match('exp*', col)]]
    data = data.apply(pd.to_numeric, errors='ignore')
    #normalized_log10_data=equalize_medians(data,log10=False)
    normalized_log10_data = data;

    ### histograms
#    plt.figure();
#    fig=plot_histograms(data)
#    print_svg()	
#    plt.close()
    if 'impute' not in plot_type:
        print("Content-Type: image/svg+xml/text/html")
        print()  # Blank line indicating the end of headers

    if 'impute' in plot_type:
        result = knn_impute(proteins).tolist()
        print(json.dumps(result))
    elif 'box_plot' in plot_type:
        plt.figure();
        fig=plot_boxplot(data,log10=False)
        plt.title('protein intensity (1og2)')
        print_svg('svg_box_plot')
        plt.close()
    elif 'density_plot' in plot_type:
        plt.figure();
        fig=plot_densityplot(data, log10=False)
        print_svg('svn_density_plot')
        plt.close()
    elif 'box_plot_norm' in plot_type:
        plt.figure();
        fig=plot_boxplot(normalized_log10_data, log10=False)
        plt.title('normalized protein intensity (log2)')
        print_svg('svg_norm_box_plot')
        plt.close()
    elif 'cor_heatmap' in plot_type:
        plt.figure()
        cross_corr2 = plot_sample_correlations_sns(normalized_log10_data,
						data_columns="exp.*",
						mode="heatmap",
						log10=False)
        print_svg('svg_cor_heatmap')
        plt.close()
#        plt.figure()
#        normalized_log10_data=equalize_medians(data,log10=False)
#        cross_corr2 = plot_sample_correlations(normalized_log10_data,                                                                                                 data_columns="I_.*",
#                                                mode="heatmap",
#                                                log10=False)
#        print_svg_plotly(cross_corr2, 'svg_cor_heatmap')
#        plt.close()
#
    elif 'volc' in plot_type:
#        c1 = 'exp1'
#        c2 = 'exp2'
#        result_square, fig_square = run_ttest(normalized_log10_data, 
#                                            cols_ann=['protein'],
#                                            c1_str=c1, 
#                                            c2_str=c2,
#                                            min_fc=1, 
#                                            s0=None, 
#                                            fdr=0.01, 
#                                            n_perm=10, 
#                                            plot_fdr_line=False)

        fig_square = plot_fc(proteins, min_fc=1, fdr=0.01)
        #svg_bytes = pio.to_image(fig_square, format='svg')
        #print (fig_square)
        print_svg_plotly(fig_square)

        #result_volc=sort_dataframe_result(result_volc,'fc') 
        #with pd.ExcelWriter('result_volc.xlsx', engine='xlsxwriter') as writer:
            #result_volc.to_excel(writer, index=False, sheet_name='Sheet1')
         
    else:
        print (f"ERROR: unsupported plot format {plot_type}")
    return

#### For command line usage
if __name__ == "__main__":
    main()


