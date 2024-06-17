#!/usr/local/bin/perl 
use strict;
use DBI;
use Getopt::Long;
use File::Basename;
use Cwd qw( abs_path );
use FindBin;
use Data::Dumper;
use FileHandle;

use lib "$ENV{SBEAMS}/lib/perl";

use SBEAMS::Connection;
use SBEAMS::Connection::Settings;
use SBEAMS::Connection::Tables;
use SBEAMS::PeptideAtlas;
use SBEAMS::PeptideAtlas::Settings;
use SBEAMS::PeptideAtlas::Tables;


$|++; # don't buffer output

my $sbeams = SBEAMS::Connection->new();
$sbeams->Authenticate();

my $atlas = SBEAMS::PeptideAtlas->new();
$atlas->setSBEAMS( $sbeams );

my $args = process_args();

my $atlas_build_id = $args->{atlas_build_id};
my $loc = $args->{file_location};
my $sample_id = $args->{sample_id} || '';
my $test =  $args->{test};
my $verbose = $args->{verbose};
my $msg = $sbeams->update_PA_table_variables($atlas_build_id);

# TODO
{ # Main 

  opendir (DIR, $loc);
  my @dir = readdir DIR;
  foreach my $file (@dir) {
		next if (($file eq '.') || ($file eq '..'));
		if ($file !~ /sdrf/i){
       print "skip $file\n";
       next;
    }
    open (IN, "<$loc/$file") or die "cannot open $loc/$file\n";
    my $header = <IN>;
    my %column_names =();
    my %selected_columns=(); 
    my %run2sample =();
    my %data =('run_name' => '',
               'biological_replicate' => '',
               'technical_replicate' => '',
               'fraction' => '',
               'label' => '',
               'enrichment'=> '',
               'treatment' => '',
               'growth_condition'=> '',
               'disease' =>'',
               'growth_stage' => '');
    my $pxd_id = '';
    my $label = '';
    if ($file =~ /(PXD\d+)/){
      $pxd_id = $1;
      $label = $pxd_id;
      %run2sample = get_run_sample_id(pxd=>$pxd_id);
    }else{
      if (! $sample_id){
        die "no PXD identifier in inputfile, need to provide sample_id\n";
      }else{
        $label=$sample_id;
        %run2sample = get_run_sample_id(sample_id=>$sample_id);
      }
    }
    print "\nsample_id=$sample_id pxd=$pxd_id\n";

    $header=~s/[\r\n]//g;
    my $idx=0;
    foreach my $col (split(/\t/, $header)){
      if ($col =~ /(data file|biological replicate|technical replicate|fraction identifier|label|enrichment|treatment|growth condition|growth stage|disease)/i){
					if ($col =~ /^.*\[(.*)\]$/){
            $col = $1;
          }
          if ($col eq 'data file'){
            $col = 'run_name';
          }elsif($col =~ /enrichment.*/){
            $col = 'enrichment';
          }elsif ($col =~ /fraction.*/i){
            $col = 'fraction';
          }
          $col =~ s/^\s+//;
          $col =~ s/\s+$//;
          $col =~ s/\s+/_/g;
          ## use first value in table
          if (not defined $column_names{$col}){
            $column_names{lc($col)} = $idx;
            $selected_columns{$idx} = lc($col);
          }
      }
      $idx++;
    }
    LOOP:foreach my $line (<IN>){
      $line =~s/[\r\n]//g;
      my $idx=0;
      foreach my $col (split(/\t/, $line)){
        if ($selected_columns{$idx}){
           if ($selected_columns{$idx} eq 'run_name'){
             $col =~ s/^(.*)\..*/$1/;
             $col =~ s/[\#\s+]/_/g;
             $col =~ s/.wiff.*//;
           }
           ## make sure have int value in biological_replicate, technical_replicate, fraction column
           if ($selected_columns{$idx} =~ /(biological_replicate|technical_replicate|fraction)/){
             if ($col =~ /(Not available|none|NA)/i){
               $col = '';
             }elsif ($col =~ /\D/){
								print "$label $line\tbiological_replicate, technical_replicate, fraction need to have integer value\n";
								next LOOP;
						 }
           }
           my $column_name = $selected_columns{$idx};
           if (defined $data{$column_name}){
             $data{$column_name} = $col;
           }else{
             print "$label sample_metadata table dont't have $column_name column\n";
           }
        }
        $idx++;  
      }
			if ($data{'run_name'}){
				my $sample_id = $run2sample{$label}{$data{run_name}};
				$data{sample_id} = $sample_id;
				if (! $sample_id){
					 print "WARNING: sample_id not found for $data{run_name}, skip\n";
					 next;
				}
				my $sql = qq~
						 SELECT SAMPLE_METADATA_ID 
						 FROM  $TBAT_SAMPLE_METADATA
						 WHERE SAMPLE_ID = $sample_id 
						 and run_name = '$data{run_name}'
				~;
				my @ids = $sbeams->selectOneColumn($sql);
				if (@ids > 1){
					print "ERROR: more than one sample_metadata_id found for run=$data{run_name} sample_id=$sample_id\n";
				}elsif (! @ids){
						my $id = $sbeams->updateOrInsertRow( insert => 1,
																		table_name  => $TBAT_SAMPLE_METADATA,
																		rowdata_ref => \%data,
																		verbose     => $verbose,
																	 PK => 'SAMPLE_METADATA_ID',
																		testonly    => $test);
				}else{  

						my $success = $sbeams->updateOrInsertRow( update => 1,
																		table_name  => $TBAT_SAMPLE_METADATA,
																		rowdata_ref => \%data,
																		verbose     => $verbose,
																		 PK => 'SAMPLE_METADATA_ID',
																		PK_value    => $ids[0],
																		testonly    => $test);
				}
			}
   }
  }
} # End Main



##################################################
sub get_run_sample_id{
  my %args = @_;
  my $pxd = $args{pxd} || '';
  my $sample_id = $args{sample_id} || '';
  my %data = ();
  my $label = '';
  my $sql = qq~
    SELECT S.spectrum_name, SMP.sample_id
    FROM $TBAT_ATLAS_BUILD AB
    JOIN $TBAT_ATLAS_BUILD_SEARCH_BATCH ABSB ON (AB.ATLAS_BUILD_ID =ABSB.ATLAS_BUILD_ID)
    JOIN $TBAT_SAMPLE SMP ON ( SMP.SAMPLE_ID = ABSB.SAMPLE_ID )
    JOIN $TBAT_SPECTRUM S ON ( S.SAMPLE_ID = SMP.SAMPLE_ID )
    WHERE AB.atlas_build_id = $atlas_build_id
  ~;

  if ($pxd){
    $label = $pxd;
    $sql .= "AND SMP.REPOSITORY_IDENTIFIERS = '$pxd'";
  }elsif ($sample_id){
    $sql .= "AND SMP.sample_id = $sample_id";
    $label= $sample_id;
  } 
  my @results = $sbeams->selectSeveralColumns($sql);
  foreach my $row(@results){
    my ($spectrum_name, $id) = @$row;
    $spectrum_name =~ /^(.*)\.\d+\.\d+\.\d+$/;
    my $run_name = $1;
    $data{$label}{$run_name}= $id;
    #print "$label $run_name $id\n";
  }
  return %data;

}


sub process_args {

  my %args;
  GetOptions( \%args, 'file_location:s', 'atlas_build_id:i','sample_id:i', 'test', 'verbose' );

  print_usage() if $args{help};

  my $missing = '';
  for my $opt ( qw(file_location atlas_build_id ) ) {
    $missing = ( $missing ) ? $missing . ", $opt" : "Missing required arg(s) $opt" if !$args{$opt};
  }
  print_usage( $missing ) if $missing;
  return \%args;
}


sub print_usage {
  my $msg = shift || '';
  my $exe = basename( $0 );

  print <<"  END";
      $msg

usage: $exe \\
       -i 576 \\
       -l /proteomics/peptideatlas/archive/build_outputs/Ecoli_2024-03/sdrf   
   -h, --help           print this usage and exit
   -a, --atlas_build_id
	 -f, --file_location 
   -s, --sample_id 
   --test 
   --verbose
  END

  exit;

}

