#!/usr/local/bin/perl

###############################################################################
# Program     : GetResultset.cgi
# Author      : Eric Deutsch <edeutsch@systemsbiology.org>
# $Id: GetResultSet.cgi 6656 2011-01-24 19:04:50Z dcampbel $
#
# Description : This CGI program dumps the ResultSet data to the user
#               in various formats
#
# SBEAMS is Copyright (C) 2000-2014 Institute for Systems Biology
# This program is governed by the terms of the GNU General Public License (GPL)
# version 2 as published by the Free Software Foundation.  It is provided
# WITHOUT ANY WARRANTY.  See the full description of GPL terms in the
# LICENSE file distributed with this software.
###############################################################################


###############################################################################
# Set up all needed modules and objects
###############################################################################
use strict;
use Getopt::Long;
use FindBin;
use Data::Dumper;
use JSON;
use Encode qw(encode decode);

use lib "$FindBin::Bin/../../lib/perl";
use vars qw ($sbeams $sbeamsMOD $q $current_contact_id $current_username
             $PROG_NAME $USAGE %OPTIONS $QUIET $VERBOSE $DEBUG $TESTONLY
             $TABLE_NAME $PROGRAM_FILE_NAME $CATEGORY $DB_TABLE_NAME
             @MENU_OPTIONS);

use SBEAMS::Connection qw($q $log);
use SBEAMS::Connection::Settings;
use SBEAMS::Connection::Tables;
use SBEAMS::PeptideAtlas;

$sbeams = new SBEAMS::Connection;
my $atlas = new SBEAMS::PeptideAtlas;
$atlas->setSBEAMS( $sbeams );

my $pid = $$;


###############################################################################
# Set program name and usage banner for command line use
###############################################################################
$PROG_NAME = $FindBin::Script;
$USAGE = <<EOU;
Usage: $PROG_NAME [OPTIONS] key=value key=value ...
Options:
  --verbose n         Set verbosity level.  default is 0
  --quiet             Set flag to print nothing at all except errors
  --debug n           Set debug flag to level n
  --testonly          Set testonly flag which simulates INSERTs/UPDATEs only

 e.g.:  $PROG_NAME --verbose 2 keyword=value

EOU

#### Process options
unless (GetOptions(\%OPTIONS,"verbose:s","quiet","debug:s","quiet")) {
  print "$USAGE";
  exit;
}

$VERBOSE = $OPTIONS{"verbose"} || 0;
$QUIET = $OPTIONS{"quiet"} || 0;
$DEBUG = $OPTIONS{"debug"} || 0;
$TESTONLY = $OPTIONS{"testonly"} || 0;
if ($DEBUG) {
  print "Options settings:\n";
  print "   VERBOSE = $VERBOSE\n";
  print "     QUIET = $QUIET\n";
  print "     DEBUG = $DEBUG\n";
  print "  TESTONLY = $TESTONLY\n";
}


###############################################################################
# Set Global Variables and execute main()
###############################################################################
main();
exit(0);


###############################################################################
# Main Program:
#
# Call $sbeams->Authenticate() and exit if it fails or continue if it works.
###############################################################################
sub main {

  #### Do the SBEAMS authentication and exit if a username is not returned
  exit unless ($current_username = $sbeams->Authenticate(
    #permitted_work_groups_ref=>['xxx','yyy'],
    #connect_read_only=>1,
    allow_anonymous_access=>1,
  ));

  my $mem = $sbeams->memusage( pid => $pid );
  $log->debug( "Init: " . $mem );

  #### Read in the default input parameters
  my %parameters;
  my $n_params_found = $sbeams->parse_input_parameters(
    q=>$q,parameters_ref=>\%parameters);

  #### Decide what action to take based on information so far
  if (defined($parameters{action}) && $parameters{action} eq "???") {
    # Some action
  } else {
    #$sbeams->display_page_header();
    handle_request(ref_parameters=>\%parameters);
    #$sbeams->display_page_footer();
  }
} # end main



###############################################################################
# Handle Request
###############################################################################
sub handle_request {
  my %args = @_;

  #https://db.systemsbiology.net/devZS/sbeams/cgi/PeptideAtlas/read_query_resultset?set_name=query_guest_20240612-091802-277
  #https://db.systemsbiology.net/devZS/sbeams/cgi/PeptideAtlas/ReadQueryResult.cgi?rs_page_number=2&rs_page_size=50&rs_set_name=query_guest_20240612-091802-277&TABLE_id=TBL_0
  #
  #### Process the arguments list
  my $ref_parameters = $args{'ref_parameters'}
    || die "ref_parameters not passed";
  my %parameters = %{$ref_parameters};

  #### Define some general variables
  my ($i,$element,$key,$value,$line,$result,$sql);
  my %resultset = ();
  my $resultset_ref = \%resultset;

#  my $msg = "empirical CE is $parameters{empirical_ce}<BR>\n";
  my $rs_page_number = $parameters{rs_page} || 1; 
  my $rs_page_size = $parameters{rs_page_size} || 50;

  #### verify that needed parameters were passed
  unless ($parameters{rs_set_name}) {
    print "ERROR: not all needed parameters were passed.  This should never ".
      "happen!  Please report this error.<BR>\n";
    return;
  }

  my $rs_set_name = $parameters{rs_set_name} ;
  my $json = new JSON; 
	$json = $json->pretty([1]);

  #### Read in the result set
  $sbeams->readResultSet(  resultset_file => $parameters{rs_set_name},
                            resultset_ref => $resultset_ref,
                     query_parameters_ref => \%parameters,
                                      debug => 0 );

  my %result = ();

	#### Default format is Tab Separated Value
	my $content_type = "Content-type: text/tab-separated-values\n\n";

													# will need tweaking if we add xml.
	print $content_type;


  if ($parameters{action} eq 'total'){
     $result{totalItems} = scalar @{$resultset_ref->{data_ref}};
  }else{
		my $column_titles_ref = $parameters{'__column_titles'};
		my %column_title_hash = ();
    my $url_cols_ref = {};
		my $idx =0;

		foreach my $column_name(@$column_titles_ref){
			$column_title_hash{$column_name} =$idx;
			$idx++;
		} 
		if ($parameters{url_cols_ref}){
			$url_cols_ref = $parameters{url_cols_ref};
		}

		my $start_idx = ($rs_page_number -1) * $rs_page_size;
		my $end_idx = $start_idx + $rs_page_size;
		my @data =();

		for(my $i=0; $i< scalar @{$resultset_ref->{data_ref}}; $i++){
		  if ($i >= $start_idx &&  $i<$end_idx){
			 	if (%column_title_hash){
					 my @values = (); 
					 ## order values according to column idx;
					 foreach my $col(@$column_titles_ref){
              my $val = $resultset_ref->{data_ref}->[$i]->[$column_title_hash{$col}];
              if ($url_cols_ref->{$col}){
                 my $url = "https://db.systemsbiology.net$url_cols_ref->{$col}";
                 $url =~ s/%0V/$val/;
                 $val = "<a href='$url'>$val</a>";
              }    
					  	push @values, $val; 
					 }
					 push @data, \@values;
				}else{
					 push @data, $resultset_ref->{data_ref}[$i];
				}
		  }
       
		}
    $result{data} = \@data; 
  }
	print  $json->encode($result{data});

} # end handle_request


#sub get_tsv_format {
#
#  my $resultset_ref = shift;
#
#  my $nrows = scalar(@{$resultset_ref->{data_ref}});
#  my @tsv = ( $resultset_ref->{column_list_ref} );
#
#	for ( my $i = 0; $i < $nrows; $i++ ) {
#		my @row = @{$resultset_ref->{data_ref}->[$i]};
#		my @return_row;
#		for my $item ( @row ) {
#			$item =~ s/<[^>]*>//gm;
#			$item =~ s/\&nbsp\;//gm;
#			push @return_row, $item;
#		}
#		push @tsv, \@return_row;
#	}
#  return \@tsv;
#}
