#!/usr/local/bin/perl 

use strict;
use Getopt::Long;
use lib( "/net/db/projects/PeptideAtlas/pipeline/bin/lib" );
use FAlite;
use FindBin;
$|++;

use lib "$ENV{SBEAMS}/lib/perl/";
use vars qw ($sbeams $sbeamsMOD $current_contact_id $current_username
             $PROG_NAME $USAGE %OPTIONS $QUIET $VERBOSE $DEBUG $TESTONLY 
             $TABLE_NAME ); 

use SBEAMS::Connection;
use SBEAMS::Connection::Settings;
use SBEAMS::Connection::Tables;

$sbeams = new SBEAMS::Connection;


use SBEAMS::PeptideAtlas;
use SBEAMS::PeptideAtlas::Settings;
use SBEAMS::PeptideAtlas::Tables;

$sbeamsMOD = new SBEAMS::PeptideAtlas;
$sbeamsMOD->setSBEAMS($sbeams);

$PROG_NAME = $FindBin::Script;
$USAGE = <<EOU;
Usage: $PROG_NAME [OPTIONS]
Options:
  --verbose n                 Set verbosity level.  default is 0
  --quiet                     Set flag to print nothing at all except errors
  --debug n                   Set debug flag
  --testonly                  If set, rows in the database are not changed or added
  --load 
  --update
  --purge
  --atlas_build_id            required
  --infile                    required 
  --outfile                    output file name 
 e.g.:  ./$PROG_NAME --atlas_build_id 347 --load --infile protein_PTM_summary.txt  
        ./$PROG_NAME --atlas_build_id 437 --infile protein_PTM_summary.txt --update
        ./$PROG_NAME --atlas_build_id 437 --infile protein_PTM_summary.txt --outfile ptm_summary.txt
EOU

#### Process options
unless (GetOptions(\%OPTIONS,"verbose:s","quiet","debug:s","testonly",
        "update","atlas_build_id:i", "load", "purge", "infile|i:s", "outfile|o:s"
    )) {

    die "\n$USAGE";

}


$VERBOSE = $OPTIONS{"verbose"} || 0;
$QUIET = $OPTIONS{"quiet"} || 0;
$DEBUG = $OPTIONS{"debug"} || 0;
$TESTONLY = $OPTIONS{"testonly"} || 0;
my $LOAD = $OPTIONS{"load"} || 0;
my $PURGE = $OPTIONS{"purge"} || 0;
my $UPDATE = $OPTIONS{"update"} || 0;
my $file = $OPTIONS{"infile"} || die $USAGE;
my $outfile = $OPTIONS{"outfile"} || '';
my $atlas_build_id = $OPTIONS{"atlas_build_id"} || die $USAGE;

$sbeams->update_PA_table_variables ($atlas_build_id);

main();
exit;



###############################################################################
## usage
################################################################################

sub main {
 
  if ($LOAD || $PURGE){
     purge_table (atlas_build_id => $atlas_build_id);
  }
  my %biosequence_ids = ();

  if ($LOAD || $UPDATE || $outfile){
     get_biosequence_ids (atlas_build_id => $atlas_build_id,
                          biosequence_ids => \%biosequence_ids);
     load_table ( file => $file,
                  outfile => $outfile,
                  atlas_build_id => $atlas_build_id,
                  biosequence_ids => \%biosequence_ids,
                  update => $UPDATE);
   }
}

sub purge_table {
  my %args  = @_;
  my $atlas_build_id = $args{atlas_build_id};
  #my $sql = "DELETE FROM $TBAT_PTM_SUMMARY WHERE ATLAS_BUILD_ID = $atlas_build_id"; 
  #print "$sql\n";
  #$sbeams->do( $sql);
   my $sql = qq~
    DECLARE \@count int
    SET \@count = 100000
    WHILE \@count > 0
    BEGIN
    DELETE
    FROM $TBAT_PTM_SUMMARY
    WHERE  ID in (
      SELECT TOP (\@count) ID
      from $TBAT_PTM_SUMMARY
      WHERE atlas_build_id = $atlas_build_id
    )
    SET \@count=\@\@ROWCOUNT;
    END
   ~;

   print "Purging BIOSEQUENCE_ID_ATLAS_BUILD_SEARCH_BATCH table ...\n";
   $sbeams->executeSQL($sql);

} 
sub get_biosequence_ids{
  my %args  = @_;
  my $atlas_build_id = $args{atlas_build_id};
  my $biosequence_ids = $args{biosequence_ids};
  my $sql = qq~
    SELECT BIOSEQUENCE_NAME, BIOSEQUENCE_ID 
    FROM $TBAT_BIOSEQUENCE BS
    JOIN $TBAT_ATLAS_BUILD AB ON (BS.BIOSEQUENCE_SET_ID = AB.BIOSEQUENCE_SET_ID )
    WHERE AB.ATLAS_BUILD_ID = $atlas_build_id 
  ~;
  %$biosequence_ids = $sbeams->selectTwoColumnHash($sql);
 
}
sub load_table{  
  my %args  = @_;
  my $atlas_build_id = $args{atlas_build_id} || die "need atlas_build_id\n";
  my $biosequence_ids = $args{biosequence_ids} || die "need biosequence_ids\n"; 
  my $update = $args{update};
  my $file = $args{file} || die "need file\n";
  my $outfile = $args{outfile} || '';

  #purge_table($atlas_build_id);

  open (IN,"<$file") or die "cannot open $file\n";
  my $line = <IN>;
  my $counter = 1;
  my %results =();
  if (! $outfile){
		my $sql  = qq~
			SELECT biosequence_id, offset, ID 
			FROM $TBAT_PTM_SUMMARY 
			WHERE ATLAS_BUILD_ID = $atlas_build_id
		~;
		my @rows = $sbeams->selectSeveralColumns($sql);
		foreach my $row (@rows){
			my ($bid,$offset,$id) = @$row;
			$results{$bid}{$offset} = $id;
		}
  }else{ 
    open (OUT, ">$outfile") or die "cannot open $outfile\n";
  } 
  while ($line = <IN>){
    chomp $line;
    my @cols = split("\t", $line, -1);
    my $protein = $cols[0];
    my $biosequence_id;
    my $str = join (",", @cols[2..$#cols]); 
    #next if ($str !~ /[1-9]/);
    if (defined $biosequence_ids->{$protein}){
      $biosequence_id = $biosequence_ids->{$protein};
    }else{
      print "WARNING: cannot find biosequence id for $protein\n";
      next;
    }
    
    if ($cols[3] > 0 && ! $cols[20]){
      die "$protein $cols[1] peptide empty\n";
    }
   # if (@cols == 21){
   #  $cols[21] = 'STY:79.9663';
   # }
    if ($outfile){
      print OUT "$counter\t$atlas_build_id\t$biosequence_id\t$cols[1]\t$cols[3]\t$cols[4]\t$cols[5]\t". 
                 join("\t", @cols[7..13])."\t$cols[18]\t$cols[19]\t$cols[2]\t$cols[6]\t$cols[20]\t".
                 "$cols[14]\t$cols[16]\t$cols[17]\t$cols[15]\t$cols[21]$cols[22]\n";
      $counter++;
      print "$counter..." if ($counter %1000 == 0);
      next;
    }
    my %rowdata = (
        atlas_build_id => $atlas_build_id,
				biosequence_id => $biosequence_id,
				offset  => $cols[1],
				residue  => $cols[2],
				nObs     => $cols[3],
				one_site  => $cols[4],
				two_sites  => $cols[5],
				over_two_sites => $cols[6],
				nP01 => $cols[7],
				nP05 => $cols[8],
				nP19 => $cols[9],
				nP81 => $cols[10],
				nP95 => $cols[11],
				nP99 => $cols[12],
				nP100    => $cols[13],
        noChoice =>$cols[14],
        enrichedWithMod => $cols[15],
        enrichedNonMod => $cols[16],
        nonEnriched => $cols[17],
				isInUniProt   => $cols[18],
				isInNeXtProt => $cols[19],
        most_observed_ptm_peptide => $cols[20],
        ptm_type => $cols[21],
        category => $cols[22]
    );
    my $offset = $cols[1];
    my $success;
    if ($update){
      if (defined $results{$biosequence_id}{$offset}){ 
        $success = $sbeams->updateOrInsertRow(
            update=>1,
            table_name=>$TBAT_PTM_SUMMARY,
            rowdata_ref=>\%rowdata,
            PK => 'id',
            PK_value => $results{$biosequence_id}{$offset}, 
            verbose=>$VERBOSE,
            testonly=>$TESTONLY,
        );
     }else{
       $success = $sbeams->updateOrInsertRow(
						table_name=>$TBAT_PTM_SUMMARY,
						insert=>1,
						rowdata_ref=>\%rowdata,
						PK => 'id',
						return_PK=>1,
						verbose=>$VERBOSE,
						testonly=>$TESTONLY,
				);
     }
   }else{
       $success = $sbeams->updateOrInsertRow(
            table_name=>$TBAT_PTM_SUMMARY,
            insert=>1,
            rowdata_ref=>\%rowdata,
            PK => 'id',
            return_PK=>1,
            verbose=>$VERBOSE,
            testonly=>$TESTONLY,
        );
   }
   $counter++;
   print "$counter..." if ($counter %1000 ==0);
  }
  close OUT;
}
