package SBEAMS::PeptideAtlas::Statistics;

use strict;
use warnings;
use Data::Dumper;
use SBEAMS::Connection qw($q $log);
use SBEAMS::Connection::Tables;
use SBEAMS::Connection::Settings;
use CGI qw(:standard);

my $sbeams;


use lib "/net/dblocal/www/html/devZS/sbeams/lib/perl/SBEAMS/Extern/lib/perl5";
use Statistics::TTest;
use Statistics::Multtest qw(:all);
use Statistics::Distributions qw(tprob);
use List::Vectorize;

###############################################################################
# Constructor
###############################################################################
sub new {
  my $this = shift;
  my $class = ref($this) || $this;
  my $self = {};
  bless $self, $class;
  $sbeams = $self->getSBEAMS();
  return($self);
} # end new

###############################################################################
# getSBEAMS: Provide the main SBEAMS object
#   Copied from AtlasBuild.pm
###############################################################################
sub getSBEAMS {
  my $self = shift;
  return $sbeams || SBEAMS::Connection->new();
} # end getSBEAMS


sub ttest_analysis {
  my $self = shift;
  my %args = @_;
  my @data = @{$args{data}};
  my $results = $args{results};
  my @sample1_cols  = @{$args{sample1_cols}};
  my @sample2_cols  = @{$args{sample2_cols}};
  my $sample1_tag  = $args{sample1_tag};
  my $sample2_tag  = $args{sample2_tag};

  my @header = @{$args{column_titles_ref}};
 
	my @proteins = map { $_->[0] } @data;

	# Initialize arrays for storing results
	my (@pvalues);

	# Iterate over proteins and perform analysis
	eval {
		for my $i (0 .. $#proteins) {
				my $protein = $proteins[$i];
				my @sample1_values = ();
				my @sample2_values = (); 
				foreach my $col (@sample1_cols){
					push @sample1_values, $data[$i]->[$col] if ($data[$i]->[$col]);
				}
				foreach my $col (@sample2_cols){
					push @sample2_values, $data[$i]->[$col] if ($data[$i]->[$col]);
				}

				# Perform t-test
				my $ttest = Statistics::TTest->new();
				$ttest->load_data(\@sample1_values, \@sample2_values);
				my $tval = ( $ttest->t_statistic < 0.001)? sprintf("%.4E",$ttest->t_statistic) : sprintf("%.3f",$ttest->t_statistic); 

				#provides the probability of observing a value at least as extreme as the calculated t-statistic under the null hypothesis 
				my $pval1 = sprintf("%.4f", $ttest->t_prob);

				# Calculate degrees of freedom
				#my $df = degrees_of_freedom (\@sample1_values, \@sample2_values);
				#assume equal variances for the two samples.
				my $df = scalar(@sample1_values) + scalar(@sample2_values) - 2;

				# Calculate p-value (two-tailed)
				my $pval2 = 2*tprob($df, abs($tval));

				# Check if p-value is defined and valid
				if (defined $pval1 && $pval1 =~ /^[0-9.]+$/) {
						push @pvalues, $pval1;
				} else {
						push @pvalues, 1; # Assign a default high p-value for invalid cases
				}

				# Calculate fold change
				my $mean1 = get_mean(\@sample1_values);
				my $mean2 = get_mean(\@sample2_values);

				my $log2_fold_change = sprintf("%.3f", $mean1-$mean2);
				#print "$protein $log2_fold_change tval=$tval pval=$pval<BR>";
				# Store results
				my $input_values = "$sample1_tag(". join(",", @sample1_values) . "); $sample2_tag(" . 
													 join(",", @sample2_values) . ')';

				my @data = ($protein, "$sample1_tag-$sample2_tag", "$input_values", $log2_fold_change, $tval);
				push @data, ($pval2 < 0.001 && $pval2 > 0)? sprintf("%.3E",$pval2) : sprintf("%.3f",$pval2); 
				push @$results, \@data;
		}
  };
  if ($@) {
    return $@;
  }

  ## value are higher than other method
	# Apply Benjamini-Hochberg correction using Statistics::Multtest
	my $qvalues_ref = BH(\@pvalues);

  ## slow
  # Apply Storey and Tibshirani correction using Statistics::Multtest
  # The value of the tuning parameter 
#  my %setup = (lambda => multiply(seq(0, 50, 5), # The value of the tuning parameter 
#                         0.01), # to estimate pi_0. Must be in [0,1).
#               robust => 0); # An indicator of whether it is desired to make the estimate 
#                             ## more robust for small p-values and a direct finite sample 2stimate of pFDR
#  my $qvalues_ref =  qvalue(\@pvalues, %setup); 

	# Print results
	for my $i (0 .. scalar @$results -1) {
			my $result = $results->[$i];
      if (defined $qvalues_ref->[$i]){
        my $qval = $qvalues_ref->[$i];
        if ($qval > 0){ 
          push @{$results->[$i]},($qval < 0.001)? sprintf("%.3E",$qval) : sprintf("%.3f",$qval);
        }else{
          push @{$results->[$i]},0;
        }
      }else{
        push @{$results->[$i]}, 'NA';
      }
			#printf "%s: Fold Change: %.3f, t-value: %.3f, p-value: %.4f, q-value (FDR): %s\n",
			#		$result->{protein}, $result->{log2_fold_change}, $result->{tval}, $result->{pval}, $qval;
	}
  return 1;
}
# Function to calculate the degrees of freedom for two samples
sub degrees_of_freedom {
    my ($data1_ref, $data2_ref) = @_;
    my @data1 = @$data1_ref;
    my @data2 = @$data2_ref;

    my $variance1 = get_stddev($data1_ref);
    my $variance2 = get_stddev($data2_ref);
    my $n1 = scalar @data1;
    my $n2 = scalar @data2;

    my $numerator = (($variance1 / $n1) + ($variance2 / $n2)) ** 2;
    my $denominator = (($variance1 / $n1) ** 2 / ($n1 - 1)) + (($variance2 / $n2) ** 2 / ($n2 - 1));

    return $numerator / $denominator;
}

# Helper function to calculate mean
sub get_mean {
    my ($data) = @_;
    my $total = 0;
    my $count = 0;
    foreach my $value (@$data){
      $total += $value if (defined $value && $value ne '');
      $count++;
    }
    return $total / $count;
}

# Function to calculate the variance of an array
sub get_stddev {
    my ($data) = @_;
    my $mean = get_mean($data);
    return undef unless defined $mean;
    my $sum = 0;
    my $count = 0;
    for my $value (@$data) {
      $sum += ($value - $mean) ** 2 if defined $value;
      $count++;
    }
    return $sum / ($count - 1);  # using sample variance formula
}

# Function to calculate the median of an array
sub get_median {
    my @array = @_;
    @array = grep { defined ($_) && $_ ne '' } @array;
    @array = sort {$a <=> $b} @array;
    my $count = scalar @array;
    if ($count % 2) {
      return $array[int($count / 2)];
    } else {
      return ($array[int($count / 2) - 1] + $array[int($count / 2)]) / 2;
    }
}

# Function to calculate the interquartile range of an array
sub get_iqr {
    my @array = @_;
    @array = sort {$a <=> $b} @array;
    my $count = scalar @array;
    my $q1 = get_median(@array[0 .. int($count / 2) - 1]);
    my $q3 = get_median(@array[int(($count + 1) / 2) .. $count - 1]);
    return $q3 - $q1;
}

sub transposed_data {
  my %args  = @_;
  my $data = $args{data};
  my $num_rows = scalar @$data;
  my $num_cols = scalar @{$data->[0]};
  my @transposed_data;
  for my $col (0 .. $num_cols - 1) {
		my @column_data;
		for my $row (0 .. $num_rows - 1) {
				push @column_data, $data->[$row][$col];
		}
		push @transposed_data, \@column_data;
  }
  return @transposed_data;
}
sub iqr_normalization{
  my $self = shift;
  my %args = @_;
  my @data = @{$args{data}};
  my $normalized_data = $args{normalized_data};
  # Transpose data to work with columns
	my @transposed_data = transposed_data(data => \@data);
  # Calculate the IQR for each sample (column)
  my @iqrs = map { get_iqr(@$_) } @transposed_data;
  # Calculate the median IQR across all samples
  my $median_iqr = get_median(@iqrs);
	# Normalize each value in the data
  my $num_rows = scalar @data;
  my $num_cols = scalar @{$data[0]};

	my @normalized_data;
	for my $row (0 .. $num_rows - 1) {
			my @normalized_row;
			for my $col (0 .. $num_cols - 1) {
        if ($data[$row][$col] && $data[$row][$col] ne ''){
					my $normalized_value = $data[$row][$col] / $iqrs[$col] * $median_iqr;
					push @normalized_row, sprintf("%.4f", $normalized_value);
        }else{
          push @normalized_row, '';
        }
			}
			push @$normalized_data, \@normalized_row;
	}

}

sub equal_median_normalization{
  my $self = shift;
  my %args = @_;
  my @data = @{$args{data}};
  my $normalized_data = $args{normalized_data};
  my $num_rows = scalar @data;
  my $num_cols = scalar @{$data[0]};
  # Transpose data to work with columns
  my @transposed_data = transposed_data(data => \@data);
  # Calculate the median for each sample (column)
  my @medians = map { get_median(@$_) } @transposed_data;
  #print join(",", @medians) ."<BR>";

  # Calculate the overall median of these medians
  my $overall_median = get_median(@medians);
	# Normalize each value in the dataset to make the median of each column equal to the overall median
	for my $row (0 .. $num_rows - 1) {
			my @normalized_row;
			for my $col (0 .. $num_cols - 1) {
        if ( $data[$row][$col] && $data[$row][$col] ne ''){
					my $normalized_value = $data[$row][$col] + ($overall_median - $medians[$col]);
					push @normalized_row, sprintf("%.4f", $normalized_value);
        }else{
          push @normalized_row, '';
        }
			}
			push @$normalized_data, \@normalized_row;
	}
}
# Function to standardize the log transformed data
sub vsn {
    my @data = @_;

    # Calculate mean and standard deviation
    my $mean = get_mean(\@data);
    my $stddev = get_stddev(\@data);

    # Standardize the data
    my @standardized_data = map { defined $_ && $_ ne '' ? sprintf("%.4f", ($_ - $mean) / $stddev) : '' } @data;
    return @standardized_data;
}

sub vsn_normalization{
  my $self = shift;
  my %args = @_;
  my @data = @{$args{data}};
  my $normalized_data = $args{normalized_data};
  # Transpose data to work with columns
  my @transposed_data = transposed_data(data => \@data);

  #Apply VSN to each column
  my @vsn_data_values  = map { [ vsn(@$_) ] } @transposed_data;

  @$normalized_data = transposed_data(data => \@vsn_data_values);
}

 
1;

__DATA__

# Example data with proteins as row names and multiple replicates for each sample
my @data = (
    ["Protein", "Sample1_1", "Sample1_2", "Sample1_3", "Sample2_1", "Sample2_2", "Sample2_3"],
    ["Protein1", 1.2, 1.3, 1.4, 0.9, 0.8, 0.85],
    ["Protein2", 2.3, 2.4, 2.5, 1.8, 1.75, 1.9],
    ["Protein3", 3.4, 3.5, 3.6, 2.7, 2.6, 2.8],
    ["Protein4", 4.5, 4.6, 4.7, 3.6, 3.55, 3.65],
    ["Protein5", 5.6, 5.7, 5.8, 4.5, 4.4, 4.55]
);
