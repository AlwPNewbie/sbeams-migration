package SBEAMS::PeptideAtlas::HTMLPrinter;

###############################################################################
# Program     : SBEAMS::PeptideAtlas::HTMLPrinter
# Author      : Eric Deutsch <edeutsch@systemsbiology.org>
# $Id$
#
# Description : This is part of the SBEAMS::WebInterface module which handles
#               standardized parts of generating HTML.
#
#		This really begs to get a lot more object oriented such that
#		there are several different contexts under which the a user
#		can be in, and the header, button bar, etc. vary by context
###############################################################################
use strict;
use LWP::UserAgent;
use HTTP::Request;
use Data::Dumper;
use LWP::Simple;

use vars qw($sbeams $current_contact_id $current_username $q
             $current_work_group_id $current_work_group_name
             $current_project_id $current_project_name $current_user_context_id);
use CGI::Carp qw(fatalsToBrowser croak);
use SBEAMS::Connection qw($log $q);
use SBEAMS::Connection::DBConnector;
use SBEAMS::Connection::Settings;
use SBEAMS::Connection::TableInfo;
use SBEAMS::Connection::DataTable;
use SBEAMS::Connection::GoogleVisualization;

use SBEAMS::PeptideAtlas::Settings;
use SBEAMS::PeptideAtlas::TableInfo;
use SBEAMS::PeptideAtlas::Tables;

use SBEAMS::Proteomics::Tables;
use SBEAMS::BioLink::Tables;

###############################################################################
# printPageHeader
###############################################################################
sub printPageHeader {
  my $self = shift;
  $self->display_page_header(@_);
}

my %dataset_annotation = ();

###############################################################################
# display_page_header
###############################################################################
sub display_page_header {
  my $self = shift;

  my %args = @_;

  my $navigation_bar = $args{'navigation_bar'} || "YES";

  my $sbeams = $self->getSBEAMS();
  my $project_id = $args{'project_id'} || $sbeams->getCurrent_project_id();
  my $output_mode = $sbeams->output_mode();
  my $http_header = $sbeams->get_http_header();

  #### If the output mode is interactive text, display text header
  if ($output_mode eq 'interactive') {
    $sbeams->printTextHeader();
    return;
  }
  #### If the output mode is interactive text, display text header
  if ($output_mode =~ 'xml') {
    my $xml_header = $sbeams->get_http_header( mode => 'xml', filename => 'peptide_export.xml' );
    print STDERR $xml_header;
    print $xml_header;
    return;
  }#elsif ($output_mode =~ 'tsv'){
  #  my $text_header = $sbeams->get_http_header( mode => 'tsv', filename => 'table.tsv' );
  #  print STDERR $text_header;
  #  print $text_header;
  #  return 
  #}

  #### If the output mode is not html, then we may not want a header here
  if ($output_mode ne 'html') {

    # Caller may want header printing to be handled here
    if ( $args{force_header} ) { 

      # Add CORS headers. - Not working...?
      #$http_header = "Access-Control-Allow-Origin : *\n" . $http_header;
      #print(STDERR "ACK\n");

      # Print http header
      print $http_header if $sbeams->invocation_mode() eq 'http';
    }

    return;
  }

  #### Obtain main SBEAMS object and use its http_header
  $sbeams = $self->getSBEAMS();

  if( $sbeams->isGuestUser() || $sbeams->getCurrent_username() =~ /^reviewer/i ) {
      $self->displayGuestPageHeader( @_ );
      return;
  } elsif ( $CONFIG_SETTING{PA_USER_SKIN} ) {
    $current_username = $sbeams->getCurrent_username();
    for my $skin ( split( ",", $CONFIG_SETTING{PA_USER_SKIN} ) ) {
#      $log->debug( "skin is $skin" );
      my ( $name, $value ) = split( /::::/, $skin, -1 );
#      $log->debug( "name is $name, val is $value" );
      if ( $name eq $current_username ) {
#        $log->debug( "CUSTOM! $CONFIG_SETTING{PA_USER_SKIN}" );
        $self->displayGuestPageHeader( @_, uri => $value );
        return;
      } else {
        $log->debug( "$name does not equal $current_username" );
      }
    }
  }
  $self->displayStandardPageHeader(@_);
}


###############################################################################
# displayInternalResearcherPageHeader
###############################################################################
sub displayInternalResearcherPageHeader {
  my $self = shift;
  my %args = @_;

  my $navigation_bar = $args{'navigation_bar'} || "YES";

  my $LOGOUT_URI = "$CGI_BASE_DIR/logout.cgi";

  my $LOGOUT_LINK = qq~<A HREF="$LOGOUT_URI" class="Nav_link">LOGOUT</A>~;


  #### Obtain main SBEAMS object and use its http_header
  my $sbeams = $self->getSBEAMS();
  my $http_header = $sbeams->get_http_header();
  use LWP::UserAgent;
  use HTTP::Request;
  my $ua = LWP::UserAgent->new();
  my $skinLink = 'https://peptideatlas.org';
  my $response = $ua->request( HTTP::Request->new( GET => "$skinLink/.index.dbbrowse2021.php" ) );
  my @page = split( "\n", $response->content() );
  my $skin = '';
  my $cnt=0;
  for ( @page ) {

    $cnt++;
    $_ =~ s/\<\!-- LOGIN_LINK --\>/$LOGOUT_LINK/;
    last if $_ =~ /--- Main Page Content ---/;
    $skin .= $_;
  }

  $skin =~ s/images\//sbeams\/images\//gm;
#  $skin =~ s/\/images\//\/sbeams\/images\//gm;
  $self->{'_external_footer'}=join("\n", @page[$cnt..$#page]);
  $self->{'_external_footer'} =~ s/images\//sbeams\/images\//gm;
 
  print "$http_header\n\n";
  print <<"  END_PAGE";
  $skin
  END_PAGE

  $self->printJavascriptFunctions();
  }


###############################################################################
# displayGuestPageHeader
###############################################################################
sub displayGuestPageHeader {
  my $self = shift;
  my %args = @_;

  my $navigation_bar = $args{'navigation_bar'} || "YES";

  my $LOGIN_URI = "$SERVER_BASE_DIR$ENV{REQUEST_URI}";

  if ($LOGIN_URI =~ /\?/) {
    $LOGIN_URI .= ";force_login=yes";
  } else {
    $LOGIN_URI .= "?force_login=yes";
  }
  my $LOGIN_LINK = qq~<A HREF="$LOGIN_URI" class="Nav_link">LOGIN</A>~;

  my $doctype = '<!DOCTYPE HTML>' . "\n";
  $doctype = '' unless $args{show_doctype}; 


  my $sbeams = $self->getSBEAMS();
  my $message = $sbeams->get_page_message( q => $q );
  $current_username ||= $sbeams->getCurrent_username();


  my $cswitcher = '';
  if ( -e "$PHYSICAL_BASE_DIR/lib/perl/SBEAMS/PeptideAtlas/ContextWidget.pm" ) {
    require SBEAMS::PeptideAtlas::ContextWidget;
    my $cwidget = SBEAMS::PeptideAtlas::ContextWidget->new();
    $cswitcher = $cwidget->getContextSwitcher( username => $current_username,
					       cookie_path => $HTML_BASE_DIR );
  } else {
    $log->debug(  "$PHYSICAL_BASE_DIR/lib/perl/SBEAMS/PeptideAtlas/ContextWidget.pm doesn't exist" ) 
  }


  # Use http_header from main SBEAMS object
  my $http_header = $sbeams->get_http_header();

  my $js = "<script language='javascript' src='/sbeams/usr/javascript/sorttable.js'></script>";


  my $ua = LWP::UserAgent->new();
  my $skinLink = $args{uri} || 'http://www.peptideatlas.org/.index.dbbrowse2021.php';
  my $resource = $sbeams->getSessionAttribute( key => 'PA_resource' ) || '';
  if ( $resource eq 'SRMAtlas' ) {
    $skinLink = 'http://www.srmatlas.org/.index.dbbrowse.php';
    #$skinLink = 'http://www.srmatlas.org/.index.dbbrowse2021.php';
  } elsif ( $resource eq 'DIAAtlas' ) {
    $skinLink = 'http://www.swathatlas.org/.index.dbbrowse-swa2021.php';
  }
  my $response = $ua->request( HTTP::Request->new( GET => "$skinLink" ) );
  my @page = split( "\n", $response->content() );
  if ( $args{show_doctype} ) {
    my $first = shift @page;
    if ( $first !~ /doctype/i ) {
      unshift @page, $first;
    }
    unshift @page, $doctype;
  }

  my $skin = '';
  my $cnt=0;
  my $init = ( $args{init_tooltip} ) ? $self->init_pa_tooltip() : '';
  my $css_info = $sbeams->printStyleSheet( module_only => 1 );
  my $loadscript = "$args{onload};" if $args{onload};
  $loadscript .= 'sortables_init();' if $args{sort_tables};

  $LOGIN_LINK .= "<BR><BR><BR>\n$cswitcher<BR>\n";
  for ( @page ) {
    $cnt++;

    # Login link originates in Peptide/SRM/SWATH Atlas. This section triages link
    # to work for current dev instance
    if ( $_ =~ /force_login=yes/ ) {
      my $url = $q->self_url();
      my $sep = ( $url =~ /\?/ ) ? ';' : '?';
      $url = $url . $sep . 'force_login=yes' unless $url =~ /force_login/;
      my $site_url = $_;
      $site_url =~ s/^(.*HREF=")[^"]+(".*$)/$1$url$2/ig;
      $skin .= $site_url;
      next;
    }
#   Turned of this mechanism - was working only for Peptide Atlas, and
#   this was printing a second login link.
#    $_ =~ s/\<\!-- LOGIN_LINK --\>/$LOGIN_LINK/;

    $_ =~ s/(\<[^>]*body[^>]*\>)/$1$init$css_info/;
    if($loadscript){
      $_ =~ s/<body/$js\n<body OnLoad="$loadscript self.focus();"/;
    }
    $_ =~ s/width="680"/width="100%"/;  # resultsets are often wide...
    $_ =~ s/width="550"//;              # and IE has trouble rendering these tables
    last if $_ =~ /--- Main Page Content ---/;
    $skin .= "$_\n";
  }
 
  $self->{'_external_footer'} = join("\n", '<!--SBEAMS_PAGE_OK-->', @page[$cnt..$#page]);
  $skin =~ s#/images/#/sbeams/images/#gm unless $skinLink =~ /2021/;
  $self->{'_external_footer'} =~ s#/images/#/sbeams/images/#gm;
  #$skin =~ s#/images/#/dev2/sbeams/images/#gm;
  
  if ( $args{tracker_type} ) {
    if ( $args{tracker_type} eq 'srm' ) {
#      $skin =~ s/<!-- Google Analytics.*<!-- End Google Analytics -->//gms;
#      die $skin;
    } elsif ( $args{tracker_type} eq 'swath' ) {
      if ( $skinLink !~ /swath/ ) {
        $skin =~ s/<!-- Google Analytics.*<!-- End Google Analytics -->//gms;
        $self->{'_external_footer'} =~ s/<!-- Google Analytics.*<!-- End Google Analytics -->//gms;
        my $track_link = $skinLink;
        $track_link =~ s/.index.*$/includes\/tracker.inc.php/;
        my $response = $ua->request( HTTP::Request->new( GET => "$track_link" ) );
        my $tracker = $response->content();
        $skin =~ s/<head>/<head>\n$tracker/gi;
      }
    }
  }

  print "$http_header\n\n";
  print "$skin";
  print "$args{header_info}\n" if $args{header_info};

  $self->printJavascriptFunctions();
  print $message;
}



###############################################################################
# displayStandardPageHeader
###############################################################################
sub displayStandardPageHeader {
  my $self = shift;
  my %args = @_;

  my $navigation_bar = $args{'navigation_bar'} || "YES";

  #### Obtain main SBEAMS object and use its http_header
  my $sbeams = $self->getSBEAMS();
  my $http_header = $sbeams->get_http_header();

  my $message = $sbeams->get_page_message( q => $q );

  my $doctype = '<!DOCTYPE HTML>' . "\n";
  $doctype = '' unless $args{show_doctype}; 

  print qq~$http_header
	$doctype<HTML><HEAD>
	<TITLE>$DBTITLE - $SBEAMS_PART</TITLE>
  ~;

  $self->printJavascriptFunctions();
  $self->printStyleSheet();

  my $loadscript = "$args{onload};" || '';
  $loadscript .= 'sortables_init();' if $args{sort_tables};
  my $js =  "<SCRIPT LANGUAGE=javascript SRC=\"/sbeams/usr/javascript/sorttable.js\"></SCRIPT>";
  print "$args{header_info}\n" if $args{header_info};

  #### Determine the Title bar background decoration
  my $header_bkg = "bgcolor=\"$BGCOLOR\"";
  $header_bkg = "background=\"$HTML_BASE_DIR//images/plaintop.jpg\"" if ($DBVERSION =~ /Primary/);

  print qq~
	<META http-equiv="X-UA-Compatible" content="chrome=IE8,IE=edge">
	<!--META HTTP-EQUIV="Expires" CONTENT="Fri, Jun 12 1981 08:20:00 GMT"-->
	<!--META HTTP-EQUIV="Pragma" CONTENT="no-cache"-->
	<!--META HTTP-EQUIV="Cache-Control" CONTENT="no-cache"-->
	</HEAD>

	<!-- Background white, links blue (unvisited), navy (visited), red (active) -->
  ~;
  if($loadscript){
     print qq~
       $js
	     <BODY BGCOLOR="#FFFFFF" TEXT="#000000" LINK="#0000FF" VLINK="#000080" ALINK="#FF0000" TOPMARGIN=0 LEFTMARGIN=0 OnLoad="$loadscript self.focus();">
     ~;
  }
  else{
      print qq~
              <BODY BGCOLOR="#FFFFFF" TEXT="#000000" LINK="#0000FF" VLINK="#000080" ALINK="#FF0000" TOPMARGIN=0 LEFTMARGIN=0 OnLoad="self.focus();">
     ~;
  }
  print $self->init_pa_tooltip() if $args{init_tooltip};

  print qq~
	<table border=0 width="100%" cellspacing=0 cellpadding=1>

	<!------- Header ------------------------------------------------>
	<a name="TOP"></a>
	<tr>
	  <td bgcolor="$BARCOLOR"><a href="http://db.systemsbiology.net/"><img height=64 width=64 border=0 alt="ISB DB" src="$HTML_BASE_DIR/images/dbsmlclear.gif"></a><a href="https://db.systemsbiology.net/sbeams/cgi/main.cgi"><img height=64 width=64 border=0 alt="SBEAMS" src="$HTML_BASE_DIR/images/sbeamssmlclear.gif"></a></td>
	  <td align="left" $header_bkg><H1>$DBTITLE - $SBEAMS_PART<BR>$DBVERSION</H1></td>
	</tr>

  ~;

  #print ">>>http_header=$http_header<BR>\n";

  if ($navigation_bar eq "YES") {

    print qq~
	<!------- Button Bar -------------------------------------------->
	<tr><td bgcolor="$BARCOLOR" align="left" valign="top">
	<table border=0 width="120" cellpadding=2 cellspacing=0>

	<tr><td><a href="$CGI_BASE_DIR/main.cgi">$DBTITLE Home</a></td></tr>
	<tr><td><a href="$CGI_BASE_DIR/$SBEAMS_PART/main.cgi">$SBEAMS_PART Home</a></td></tr>
	<tr><td><a href="$CGI_BASE_DIR/logout.cgi">Logout</a></td></tr>
	<tr><td>&nbsp;</td></tr>
	<tr><td>Browse Data:</td></tr>
	<tr><td><a href="$CGI_BASE_DIR/$SBEAMS_SUBDIR/GetPeptides"><nobr>&nbsp;&nbsp;&nbsp;Browse Peptides</nobr></a></td></tr>
	<tr><td><a href="$CGI_BASE_DIR/$SBEAMS_SUBDIR/BrowseBioSequence.cgi"><nobr>&nbsp;&nbsp;&nbsp;Browse Bioseqs</nobr></a></td></tr>
	<tr><td><a href="$CGI_BASE_DIR/$SBEAMS_SUBDIR/GetPeptide"><nobr>&nbsp;&nbsp;&nbsp;Get Peptide Summary</nobr></a></td></tr>
	<tr><td><a href="$CGI_BASE_DIR/$SBEAMS_SUBDIR/GetProtein"><nobr>&nbsp;&nbsp;&nbsp;Get Protein Summary</nobr></a></td></tr>
	<tr><td><a href="$CGI_BASE_DIR/$SBEAMS_SUBDIR/showPathways"><nobr>&nbsp;&nbsp;&nbsp;Pathway Search</nobr></a></td></tr>
	<tr><td><a href="$CGI_BASE_DIR/$SBEAMS_SUBDIR/GetTransitions"><nobr>&nbsp;&nbsp;&nbsp;SRM Transitions</nobr></a></td></tr>
	<tr><td><a href="$CGI_BASE_DIR/Glycopeptide/Glyco_prediction.cgi"><nobr>&nbsp;&nbsp;&nbsp;Search Glyco-Peptides</nobr></a></td></tr>
	<tr><td>&nbsp;</td></tr>
	<tr><td>Manage Tables:</td></tr>
	<tr><td><a href="$CGI_BASE_DIR/$SBEAMS_SUBDIR/ManageTable.cgi?TABLE_NAME=AT_biosequence_set"><nobr>&nbsp;&nbsp;&nbsp;Biosequence Sets</nobr></a></td></tr>
	<tr><td><a href="$CGI_BASE_DIR/$SBEAMS_SUBDIR/ManageTable.cgi?TABLE_NAME=AT_atlas_build"><nobr>&nbsp;&nbsp;&nbsp;Atlas Builds</nobr></a></td></tr>
	<tr><td><a href="$CGI_BASE_DIR/$SBEAMS_SUBDIR/ManageTable.cgi?TABLE_NAME=AT_default_atlas_build"><nobr>&nbsp;&nbsp;&nbsp;Default Builds</nobr></a></td></tr>
	<tr><td><a href="$CGI_BASE_DIR/$SBEAMS_SUBDIR/ManageTable.cgi?TABLE_NAME=AT_atlas_search_batch"><nobr>&nbsp;&nbsp;&nbsp;Atlas Search Batches</nobr></a></td></tr>
	<tr><td><a href="$CGI_BASE_DIR/$SBEAMS_SUBDIR/ManageTable.cgi?TABLE_NAME=AT_sample"><nobr>&nbsp;&nbsp;&nbsp;Samples</nobr></a></td></tr>
	<tr><td><a href="$CGI_BASE_DIR/$SBEAMS_SUBDIR/ManageTable.cgi?TABLE_NAME=BL_dbxref"><nobr>&nbsp;&nbsp;&nbsp;DB Xrefs</nobr></a></td></tr>
	<tr><td><a href="$CGI_BASE_DIR/$SBEAMS_SUBDIR/ManageTable.cgi?TABLE_NAME=organism"><nobr>&nbsp;&nbsp;&nbsp;Organisms</nobr></a></td></tr>
	<tr><td><a href="$CGI_BASE_DIR/$SBEAMS_SUBDIR/ManageTable.cgi?TABLE_NAME=AT_publication"><nobr>&nbsp;&nbsp;&nbsp;Publications</nobr></a></td></tr>
  <tr><td><a href="$CGI_BASE_DIR/$SBEAMS_SUBDIR/ManageTable.cgi?TABLE_NAME=AT_elution_time_type"><nobr>&nbsp;&nbsp;&nbsp;Elution Time Type</nobr></a></td></tr>
  <tr><td><a href="$CGI_BASE_DIR/$SBEAMS_SUBDIR/ManageTable.cgi?TABLE_NAME=AT_instrument_type"><nobr>&nbsp;&nbsp;&nbsp;Instrument Type</nobr></a></td></tr>
  <tr><td><a href="$CGI_BASE_DIR/$SBEAMS_SUBDIR/ManageTable.cgi?TABLE_NAME=AT_fragmentation_type"><nobr>&nbsp;&nbsp;&nbsp;Fragmentation Type</nobr></a></td></tr>

	<tr><td><a href="$CGI_BASE_DIR/$SBEAMS_SUBDIR/ManageTable.cgi?TABLE_NAME=AT_PASS_submitter"><nobr>&nbsp;&nbsp;&nbsp;PASS</nobr></a></td></tr>
  <tr><td><a href="$CGI_BASE_DIR/$SBEAMS_SUBDIR/ManageTable.cgi?TABLE_NAME=AT_public_data_repository"><nobr>&nbsp;&nbsp;&nbsp;Public data repository</nobr></a></td></tr>
	<tr><td>&nbsp;</td></tr>
	<tr><td>Annotations:</td></tr>
	<tr><td><a href="$CGI_BASE_DIR/$SBEAMS_SUBDIR/ManageTable.cgi?TABLE_NAME=AT_spectrum_annotation"><nobr>&nbsp;&nbsp;&nbsp;Spectra</nobr></a></td></tr>
	<tr><td><a href="$CGI_BASE_DIR/$SBEAMS_SUBDIR/ManageTable.cgi?TABLE_NAME=AT_peptide_annotation"><nobr>&nbsp;&nbsp;&nbsp;Peptides</nobr></a></td></tr>
	<tr><td><a href="$CGI_BASE_DIR/$SBEAMS_SUBDIR/ManageTable.cgi?TABLE_NAME=AT_modified_peptide_annotation"><nobr>&nbsp;&nbsp;&nbsp;Modified Peptides</nobr></a></td></tr>
	<tr><td><a href="$CGI_BASE_DIR/$SBEAMS_SUBDIR/ManageTable.cgi?TABLE_NAME=AT_spectrum_annotation_level"><nobr>&nbsp;&nbsp;&nbsp;Spectra Levels</nobr></a></td></tr>
	<tr><td><a href="$CGI_BASE_DIR/$SBEAMS_SUBDIR/ManageTable.cgi?TABLE_NAME=AT_transition_suitability_level"><nobr>&nbsp;&nbsp;&nbsp;Transition Levels</nobr></a></td></tr>
	</table>
	</td>

	<!-------- Main Page ------------------------------------------->
	<td valign=top>
	<table border=0 bgcolor="#ffffff" cellpadding=4>
	<tr><td>$message

    ~;
  } else {
      print qq~
	  </TABLE>
      ~;
  }

}

# 	<table border=0 width="680" bgcolor="#ffffff" cellpadding=4>


###############################################################################
# printStyleSheet
#
# Print the standard style sheet for pages.  Use a font size of 10pt if
# remote client is on Windows, else use 12pt.  This ends up making fonts
# appear the same size on Windows+IE and Linux+Netscape.  Other tweaks for
# different browsers might be appropriate.
###############################################################################
sub printStyleSheet {
    my $self = shift;

    #### Obtain main SBEAMS object and use its style sheet
    $sbeams = $self->getSBEAMS();
    $sbeams->printStyleSheet();
    print <<"    END";
    <STYLE>
     TD { background-repeat: no-repeat; }
    </STYLE>
    END

}


###############################################################################
# printJavascriptFunctions
#
# Print the standard Javascript functions that should appear at the top of
# most pages.  There probably should be some customization allowance here.
# Not sure how to design that yet.
###############################################################################
sub printJavascriptFunctions {
    my $self = shift;
    my $javascript_includes = shift;


    print qq~
	<SCRIPT LANGUAGE="JavaScript">
	<!--

	function refreshDocument() {
            //confirm( "apply_action ="+document.MainForm.apply_action.options[0].selected+"=");
            document.MainForm.apply_action_hidden.value = "REFRESH";
            document.MainForm.action.value = "REFRESH";
	    document.MainForm.submit();
	} // end refreshDocument


	function showPassed(input_field) {
            //confirm( "input_field ="+input_field+"=");
            confirm( "selected option ="+document.forms[0].slide_id.options[document.forms[0].slide_id.selectedIndex].text+"=");
	    return;
	} // end showPassed



        // -->
        </SCRIPT>
    ~;

}


###############################################################################
# printPageFooter
###############################################################################
sub printPageFooter {
  my $self = shift;
  $self->display_page_footer(@_);
}


###############################################################################
# display_page_footer
###############################################################################
sub display_page_footer {
  my $self = shift;
  my %args = @_;


  #### If the output mode is interactive text, display text header
  my $sbeams = $self->getSBEAMS();
  if ($sbeams->output_mode() eq 'interactive') {
    $sbeams->printTextHeader(%args);
    return;
  }


  #### If the output mode is not html, then we don't want a header here
  if ($sbeams->output_mode() ne 'html') {
    return;
  }


  #### Process the arguments list
  my $close_tables = $args{'close_tables'} || 'YES';
  my $display_footer = $args{'display_footer'} || 'YES';
  my $separator_bar = $args{'separator_bar'} || 'NO';

   if($self->{'_external_footer'}) {


	  # Have to fish for error file, since we are using Carp (can't pass sbeams
		# object ).
	  my $errfile = 'Error-' . getppid();
 	  my $is_error = 0;
 	  if ( $sbeams->doesSBEAMSTempFileExist( filename => $errfile ) ) {
			$sbeams->deleteSBEAMSTempFile( filename => $errfile );
			$is_error++;
		}

   	$self->{'_external_footer'} =~ s/SBEAMS_PAGE_OK/SBEAMS_PAGE_ERROR/gmi if $is_error;

   	print "$self->{'_external_footer'}\n";
     	return;
  }


  #### If closing the content tables is desired
  if ($close_tables eq 'YES') {
    print qq~
	</TD></TR></TABLE>
	</TD></TR></TABLE>
    ~;
  }


  #### If displaying a fat bar separtor is desired
  if ($separator_bar eq 'YES') {
    print "<BR><HR SIZE=5 NOSHADE><BR>\n";
  }


  #### If finishing up the page completely is desired
  if ($display_footer eq 'YES') {
    my $close_main_tables = 'NO';
    $close_main_tables = 'YES' unless ($close_tables eq 'YES');

    #### Default to the Core footer
    $sbeams->display_page_footer(display_footer=>'YES',
      close_tables=>$close_main_tables);
  }

}

#+
# @narg header_text - Text to render as HTML
# @narg header_element - XML friendly header element
# @narg anchor - text to create <A> anchor in HTML mode
# @narg bold - header bold text
# @narg list_items - ref to array of anon hashes of form key=>value for two column table
# @narg width - fixed pixel width for table
# @narg header_width - fixed pixel width for heading
# @narg key_width - fixed percentage of width for key
sub encodeFullSectionList {
  my $self = shift || die ("self not passed");
  my %args = ( width => 600,
               header_width => 900,
               bold => 1,
               key_width => 20,
               @_
             );

  # Default to BOLD
  unless ( $args{header_text} && $args{list_items} ) {
    $log->error( "Required parameters not supplied" );
    return '';
  }
  $args{header_element} ||= $args{header_text};

  my $buffer;
  my $sbeams = $self->getSBEAMS();
  if ( $sbeams->output_mode() =~ /html/i ) {
    $buffer = "<table style='min-width:$args{width};margin-left:15px;'>\n";
    $buffer .= $self->encodeSectionHeader( %args, text => $args{header_text} );
    for my $item ( @{$args{list_items}} ) {
      $buffer .= $self->encodeSectionItem( key => $item->{key}, value => $item->{value}, tr_info => 'class="hoverable"', key_width => $args{key_width} . "%" ) . "\n";
    } 
    $buffer .= "</table>\n";

  } elsif ( $sbeams->output_mode() =~ /xml/i ) {
    $buffer = "<$args{header_element}>\n";
    for my $item ( @{$args{list_items}} ) {
      my $key = $item->{key};
      $key =~ s/\s//g;
      $buffer .= "<li $key='$item->{value}'/>\n";
    } 
    $buffer .= "</$args{header_element}>\n";
  } else {
    $buffer = 'Sorry, yet to be implemented!';
  }
  return $buffer;
}


###############################################################################
# encodeSectionHeader
###############################################################################
sub encodeSectionHeader {
  my $METHOD = 'encodeSectionHeader';
  my $self = shift || die ("self not passed");
  my %args = @_;

  my $colspan = $args{colspan} || 2;
  
  # Default to BOLD
  $args{bold} = 1 if !defined $args{bold};

  my $text = $args{text} || '';
  $text = "<B>$text</B>" if $args{bold};

  my $link = $args{link} || '';

  $args{anchor} ||= $args{text};
  my $anchor = ( $args{anchor} ) ? "<A NAME='$args{anchor}'></A>" : '';

  my $mouseover = '';
  if ( $args{mouseover} ) {
    my $qmark = "<img src='$HTML_BASE_DIR/images/greyqmark.gif' />";
    $text = "<div title='$args{mouseover}'>$text&nbsp;$qmark</div>";
  }

#        <TR><TD colspan="2" background="$HTML_BASE_DIR/images/fade_orange_header_2.png" width="600">$link<font color="white">$anchor$text</font></TD></TR>
#        <TR><TD colspan="2" class=fade_header width="600">$link<font color="white">$anchor$text</font></TD></TR>

  my $buffer = qq~
        <TR><TD colspan="$colspan" style="background-repeat: no-repeat; background-image: url('$HTML_BASE_DIR/images/fade_orange_header_2.png')" width="600">$link<font color="white">$anchor$text</font></TD></TR>
~;

  if ( $args{LMTABS} ) {
    my $sbeams = $self->getSBEAMS();
    $args{no_toggle} ||= 0;
    $args{divname} ||= $sbeams->getRandomString( num_chars => 20 );
    $buffer = $sbeams->make_toggle_section( neutraltext => $text,
					    sticky => 1,
					    name => $args{divname},
					    tooltip => 'Show/Hide Section',
					    barlink => 1,
					    no_toggle => $args{no_toggle},
					    visible => 1,
					    content => $text
	);
  }

###        <TR><TD colspan="$colspan" $link style="background:#f3f1e4;color:#555;border-top:1px solid #b00;border-left:15px solid #b00;padding:0.5em">$anchor$text</TD></TR> e2dcc2


  return $buffer;
  return "<table>$buffer</table>";
}


###############################################################################
# encodeSectionItem
###############################################################################
sub encodeSectionItem {
  my $METHOD = 'encodeSectionItem';
  my $self = shift || die ("self not passed");
  my %args = @_;

  my $key = $args{key} || '';

  my $value = $args{value};
  $value = '' if !defined $value;

  my $desc_nowrap = ( $args{nowrap_description} ) ? 'NOWRAP' : '';

  my $url = $args{url} || '';
  my $kwid = ( $args{key_width} ) ? "WIDTH='$args{key_width}'" : '';
  my $vwid = ( $args{val_width} ) ? "WIDTH='$args{val_width}'" : '';
  
  my $tr = $args{tr_info} || ''; 

  $url =~ s/ /+/g;
  my $astart = '<b>';
  my $aend = '</b>';
  if ($url) {
    $astart = qq~<A HREF="$url" target="_blank">~;
    $aend = qq~</A>~;
  }

  my $buffer = '';
  if ($tr =~ /hoverable/) {  # new L&F?
# e2dcc2
    $buffer = qq~
        <TR $tr><TD class="key" NOWRAP $kwid>$key</TD><TD class="value" $desc_nowrap $vwid>$astart$value$aend</TD></TR>
~;
  }
  else {
    $buffer = qq~
        <TR $tr><TD NOWRAP bgcolor="cccccc" $kwid>$key</TD><TD $desc_nowrap $vwid>$astart$value$aend</TD></TR>
~;
  }


  return $buffer;

}


###############################################################################
# encodeSectionTable
###############################################################################
#
#
# @nparam  rows             Reference to array of row arrayrefs.  Required.
# @nparam  header           rows includes header info, default 0
# @nparam  width            Width of table
# @nparam  nocenter         Overrides centering of header info
# @nparam  align            Ref to array of positional info for columns -
#                           right, left, and center
# @nparam  rows_to_show     Number of rows to show initially, others are loaded
#                           into page but hidden by CSS and showable view JS.
# @nparam  nowrap           ref to array of col indexes on which to set nowrap
# @nparam  maxrows          Maximum number of rows the table can contain, 
#                           irrespective of the number displayed.
# @nparam  chg_bkg_idx      Index of field on which to trigger alternating
#                           colors.
# @nparam  bkg_interval     numeric interval on which to alternate colors, 
#                           superceded by chk_bkg_idx
# @nparam  set_download     Have table include download link if set to a true
#                           value.  If text value is used then it will be used
#                           as text for the link.
# @nparam  sortable         If true, make data table sortable
#
#
# 
#
sub encodeSectionTable {
  my $METHOD = 'encodeSectionTable';
  my $self = shift || die ("self not passed");
  my %args = @_;

  $args{rs_params} ||= {};

  my $pre_text = '';
  my $class_def = '';
  my $id = $args{table_id} || '';
  my $file_prefix = '';
  if ( $args{file_prefix}){
    $file_prefix = $args{file_prefix};
  }elsif($id){
    $file_prefix = $id ."_";
  }else{
    $file_prefix = 'mrm_';
  }

  if ( $args{sortable} ) {
    if ( !$self->{_included_sortable} ) {
      $pre_text = qq~
      <SCRIPT LANGUAGE=javascript SRC="$HTML_BASE_DIR/usr/javascript/sorttable.js"></SCRIPT>
      ~;
      $self->{_included_sortable}++;
    }
    $class_def = $args{class} || 'PA_sort_table';
  }

  my @table_attrs = ( 'border' => 0, );
  if ( $args{width} ) {
    push @table_attrs, 'WIDTH', $args{width};
  }
  my $tr_info = $args{tr_info} || 'NOOP=1 ';
  my $y_scroll = $args{y_scroll} || '';
  my $tab = SBEAMS::Connection::DataTable->new( @table_attrs,
						__use_thead => 1,
						__tr_info => $tr_info,
            __y_scroll => $y_scroll,
						CLASS => $class_def,
						ID => $id  );

  my $num_cols = 0;
  my $rs_link = '';
  my $rs_name = '';

  my @rs_data =  @{$args{rows}};

  if ( $args{set_download} ) {
    my $rs_headers = shift( @rs_data );
    $rs_headers = $args{rs_headings} if $args{rs_headings};
    
    $rs_name = $self->make_resultset( rs_data => \@rs_data, 
                                      headers => $rs_headers,
                                  file_prefix => $file_prefix,
                                   rs_params => $args{rs_params} );

    my $tsv_link = "<a href='$CGI_BASE_DIR/GetResultSet.cgi/$rs_name.tsv?rs_set_name=$rs_name&format=tsv;remove_markup=1' TITLE='Download table as tab-delimited text file'>TSV</a>";
    my @downloads = ( $tsv_link ); 
    if ( $args{download_links} ) {
      push @downloads, @{$args{download_links}};
    }
    my $download_links = join( ', ', @downloads );
    $rs_link = "<SPAN CLASS=info_box>Download as: $download_links</SPAN>";

    if ( $args{download_form} ) {
      my $hidden = qq~
      <INPUT TYPE=HIDDEN NAME=rs_set_name VALUE=$rs_name>
      <INPUT TYPE=HIDDEN NAME=remove_markup VALUE=1>
      <INPUT TYPE=HIDDEN NAME=tsv_output VALUE=1>
      ~;
      $rs_link = $args{download_form};
      $rs_link =~ s/HIDDEN_PLACEHOLDER/$hidden/m;
    }

  }

  return '' unless $args{rows};
  $args{header} ||= 0;

  $args{max_rows} ||= 0;
  my $sbeams = $self->getSBEAMS();
  my $prefix = $sbeams->getRandomString( num_chars => 8, 
                                          char_set => ['A'..'Z', 'a'..'z'] );
  my $first = 1;
  my $bgcolor = '#f3f1e4'; #C0D0C0';
  my $chg_idx;
  my $rcnt = ( $args{header} ) ? 1 : 0;
#
  my $msg = '';
  for my $row ( @{$args{rows}} ) {
    $num_cols = scalar( @$row ) unless $num_cols;

    ## add thousand separator to numbers 
    for (my $i=0; $i< scalar @$row; $i++){
      if ($row->[$i] =~ /^\d+$/){
         my $num = $row->[$i];
         while ($num =~ s/^(-?\d+)(\d{3}(?:,\d{3})*(?:\.\d+)*)$/$1,$2/){};
         $row->[$i] = $num;
      }elsif($row->[$i] =~ /^(\d+)\s([\(\)\d\.\%]+)$/){
         my $num = $1;
         my $pect = $2;
         while ($num =~ s/^(-?\d+)(\d{3}(?:,\d{3})*(?:\.\d+)*)$/$1,$2/){};
         $row->[$i] = "$num $pect";
      }
    }  

    $tab->addRow( $row );
    $rcnt++;

    if ( defined $args{chg_bkg_idx} ) { # alternate on index
      if ( !$chg_idx ) {
        $chg_idx = $row->[$args{chg_bkg_idx}];
      } elsif ( $chg_idx ne $row->[$args{chg_bkg_idx}] ) {
#       $bgcolor = ( $bgcolor eq '#C0D0C0' ) ? '#F5F5F5' : '#C0D0C0';
        $bgcolor = ( $bgcolor eq '#f3f1e4' ) ? '#d3d1c4' : '#f3f1e4';
        $chg_idx = $row->[$args{chg_bkg_idx}];
      }

    } elsif ( $args{bkg_interval} ) { # alternate on n_rows
				unless ( $rcnt % $args{bkg_interval} ) {
			#        $bgcolor = ( $bgcolor eq '#C0D0C0' ) ? '#F5F5F5' : '#C0D0C0';
					$bgcolor = ( $bgcolor eq '#f3f1e4' ) ? '#d3d1c4' : '#f3f1e4';
						}
    } elsif ( $args{bg_color} ) { # single solid color
      $bgcolor = $args{bg_color};
    }
    $tab->setRowAttr(  ROWS => [$tab->getRowNum()], BGCOLOR => $bgcolor );
    my $num_rows = $tab->getRowNum() - 1;

    if ( $args{rows_to_show} && $args{rows_to_show} < $num_rows ) {
      $tab->setRowAttr( ROWS => [$tab->getRowNum()], ID => $prefix . '_toggle', 
                                                     NAME => $prefix . '_toggle', 
                                                     CLASS => 'hidden' ); 
    }
    # Message regarding truncated rows.
    if ( $args{max_rows} && $args{max_rows} <= $num_rows ) {
      my $span = scalar( @$row );
      $msg = "Table truncated at $args{max_rows} rows";
      if  ( $args{set_download} ) {
        $msg .= "(won't affect download, " . scalar( @{$args{rows}} ) . ' total rows)';
      }
      if ( $args{truncate_msg} ) {
        $msg .= ". $args{truncate_msg}";
      }
      if ( !$args{truncate_msg_as_text} ) {
        $tab->addRow( [$sbeams->makeInfoText("<I>$msg</I>")] );
      }
      $tab->setColAttr( ROWS => [$tab->getRowNum()], 
                     COLSPAN => $span, 
                        COLS => [ 1 ], 
                       ALIGN => 'CENTER' ); 

      if ( $args{rows_to_show} && $args{rows_to_show} < $num_rows ) {
        $tab->setRowAttr( ROWS => [$tab->getRowNum()], ID => $prefix . '_toggle', 
                                                     NAME => $prefix . '_toggle', 
                                                     CLASS => 'hidden' ); 
      } 
      last;
    }
  }

  my $nowrap = $args{nowrap} || [];
  if ( scalar( @$nowrap ) ) {
    $tab->setColAttr( COLS => $nowrap, NOWRAP => 1 ); 
  }
  
  # How many do we have?
  my $tot = $tab->getRowNum();
  my $closelink;
  if ( $args{rows_to_show} && $args{rows_to_show} < $tot - 1 ) {
    $closelink = $self->add_tabletoggle_js(); 
    $closelink .= "\n<a href='#' onclick='toggle_em(\"$prefix\");return false;'><span id='${prefix}_text' name='${prefix}_text' >[Show more rows]</a>";
  }

#  # if no wrapping desired...
#  if ( $args{nowrap} ) {
#    $tab->setRowAttr( ROWS => [1..$tot], NOWRAP => 1 ); 
#  }

  # Set header attributes
  if ( $args{header} ) {
    if ( $args{sortable} ) {
#      $tab->setRowAttr( ROWS => [1], BGCOLOR => '#0000A0', CLASS => 'sortheader' );
      $tab->setRowAttr( ROWS => [1], BGCOLOR => '#003F72', CLASS => 'sortheader' );
    } else {
      $tab->setRowAttr( ROWS => [1], BGCOLOR => '#003F72', CLASS => 'sortheader' );
    }
    if ($args{header_sticky}){
      $tab->setRowAttr( ROWS => [1], STYLE=>'position:sticky;top:0;')
    }
    $tab->setRowAttr( ROWS => [1], ALIGN => 'CENTER' ) unless $args{nocenter};
  }

  if ( $args{align} ) {
    for ( my $i = 0; $i <= $#{$args{align}}; $i++ ) {
      $tab->setColAttr( ROWS => [2..$tot], COLS => [$i + 1], ALIGN => $args{align}->[$i] );
    }
  }
  if ( $args{has_key} ) {
      $tab->setColAttr( ROWS => [2..$tot], COLS => [1], CLASS => 'key' );
  }

#  $tab->addRow( [$closelink] );
  if ( $args{table_only} ) {
    return "$tab";
  }

  my $html = "$pre_text\n";
  my $help = $args{help_text} || '';

  if ( $html && $args{manual_widgets} ) {
    $html .= "<tr $tr_info><td nowrap ALIGN=left>$closelink</td><td align=center>$help</td><td nowrap ALIGN=right>$rs_link</td></tr>\n";
  }

  if ( 0 && ( $rs_link || $args{change_form} ) ) {
    
    if ( !$rs_link ) {
      $html .= "<TR><TD NOWRAP ALIGN=left>$args{change_form}</TD></TR>\n";
    } elsif ( !$args{change_form} ) {
      $html .= "<TR><TD NOWRAP COLSPAN=$num_cols ALIGN=right>$rs_link</TD></TR>\n";
    } else {
      $html .= "<TR><TD NOWRAP ALIGN=left>$args{change_form}</TD>\n";
      $html .= "<TD NOWRAP ALIGN=right>$rs_link</TD></TR>\n";
    }
  }

  my $colspan = $args{colspan} || 2;
  $html .= "<TR><TD NOWRAP COLSPAN=$colspan>$tab</TD></TR>";
  $html .= '</TABLE>' if $args{close_table};
   
  if ( wantarray ) {

    if ( $args{truncate_msg_as_text} ) {
      $html = "<tr $tr_info><td>$closelink</td></tr>" . $html;
      return ($html, $rs_name, $msg);
    } else {
      return ($html, $rs_name);
    }
  } else { 
    if ( $args{unified_widgets} ) {
      my $widget = "<tr $tr_info><td align=left>$closelink</td>";
      $widget .= "<td align=right nowrap=1>$rs_link</td>"; 
      $widget .= "</tr>";
      return( $widget . $html ); 
    } else {
      $html = "<tr $tr_info><td>$closelink</td></tr>" . $html;
      return $html;
    }
  }
}

###############################################################################
# display sample plot
###############################################################################
sub getSamplePlotDisplay {
  my $self = shift;
  my $sbeams = $self->getSBEAMS();
  my %args = @_;
  my $SUB_NAME = 'getSamplePlotDisplay';

  for my $arg ( qw( n_obs obs_per_million ) ) {
    unless ( defined ($args{$arg} ) ) {
      $log->error( "Missing required argument $arg" );
      return undef;
    }
  }

  $args{no_header} ||= 0;

  my $header = '';
  unless ( $args{no_header} ) {

    if ( $args{link} ) {
      $header = $self->encodeSectionHeader( text => 'Observed in Experiments:',
					  anchor => 'samples',
					    link => $args{link} , LMTABS=>1, no_toggle=>1);
    } else {
      $header = $self->encodeSectionHeader( text => 'Observed in Experiments:',
					  anchor => 'samples', LMTABS=>1, no_toggle=>1 );
    }
  }

  my $trinfo = $args{tr_info} || '';
  my $height = 50  + scalar( @{$args{n_obs}} ) * 12;

  # unable to hide legend as desired. Works on google viz playground, but not here (Alas!).
  my $GV = SBEAMS::Connection::GoogleVisualization->new();
  my $chart = $GV->setDrawBarChart(  samples => $args{n_obs},
                                     options => '', # qq~ legend: { position: 'none'}, title: "Number of Observations" ~, 
                                  data_types => [ 'string', 'number' ],
                                    headings => [ 'Sample Name', 'Number of Observations' ],
                                    chart_div => 'chart1_div',
                                  );
  my $chart_2 = $GV->setDrawBarChart(samples => $args{obs_per_million},
                                     options => '', # qq~ legend: { position: 'none'}, title: "Obs per million spectra" ~,
                                  data_types => [ 'string', 'number' ],
                                    headings => [ 'Sample Name', 'Obs per million spectra' ],
                                   chart_div => 'chart1_div',
                                      no_div => 1,

      );
  my $h_info = $GV->getHeaderInfo();
  $chart = qq~
    <script type='text/javascript'>
    function toggle_plot() { 
      if ( window['next_plot'] == undefined ) {
        window['next_plot'] = 'chart1_div';
      }
      if ( window['next_plot'] == 'chart1_div' ) {
        draw_chart1();
        window['next_plot'] = 'chart2_div';
        document.getElementById('toggle_button').innerHTML = 'Show Obs Per Million Spectra Plot';
      } else {
        draw_chart2();
        window['next_plot'] = 'chart1_div';
        document.getElementById('toggle_button').innerHTML = 'Show Total Obs Plot';
      }
    } 
  </script>
    $h_info
    <TR $trinfo>
      <TD></TD>
      <TD><button type='button' id='toggle_button' onclick=toggle_plot()>Show Total Obs Plot</button> 
          &nbsp;$chart_2
          &nbsp;$chart 
      </TD>
    </TR>
  ~;

  return ( wantarray() ) ? ($header, $chart) : $header . "\n" . $chart;
}

###############################################################################
# displaySampleMap
###############################################################################
sub getSampleMapDisplay {
  my $self = shift;
  my $sbeams = $self->getSBEAMS();
  my %args = ( peptide_field => 'peptide_accession', @_ );
  my $in = join( ", ", keys( %{$args{instance_ids}} ) );
  my $sample_category = $args{sample_category} || 0;
  my $sample_category_contraint = $args{sample_category_contraint} || '';

  return unless $in;

  
  my $header = '';
  if ( $args{link} ) {
    $header .= $self->encodeSectionHeader( text => 'Experiment peptide map:',
					   link => $args{link},
					   anchor => 'samples'
                                          );
  } else {
    $header .= $self->encodeSectionHeader( text => 'Experiment peptide map:',
					   anchor => 'samples'
                                         );
  }
  $header = '' if $args{no_header};

  my $html = '';
  my $trinfo = $args{tr_info} || '';
  my $order_by = 'ORDER BY sample_tag ASC';
  my $sample_tag = 'sample_tag';
  my $sample_category_clause = '';
  if ($sample_category){
    $order_by = 'ORDER BY sample_name ASC';
    $sample_tag = "SC.name + ':' + sample_tag  as sample_name";
  }
  if ($sample_category_contraint && $sample_category_contraint ne "''"){
     $sample_category_clause = "AND SC.name in ($sample_category_contraint)";
  }

  my $sql = qq~     
  	SELECT DISTINCT SB.atlas_search_batch_id, $sample_tag, 
		PISB.n_observations, $args{peptide_field},
    CASE when n_genome_locations = 1 THEN 1 ELSE 2 END
		FROM $TBAT_ATLAS_SEARCH_BATCH SB 
	  JOIN $TBAT_SAMPLE S ON s.sample_id = SB.sample_id
	  JOIN $TBAT_PEPTIDE_INSTANCE_SEARCH_BATCH PISB ON PISB.atlas_search_batch_id = SB.atlas_search_batch_id
	  JOIN $TBAT_PEPTIDE_INSTANCE PI ON PI.peptide_instance_id = PISB.peptide_instance_id
	  JOIN $TBAT_PEPTIDE P ON P.peptide_id = PI.peptide_id
    LEFT JOIN $TBAT_SAMPLE_CATEGORY SC ON (SC.id = S.sample_category_id)
    WHERE PI.peptide_instance_id IN ( $in )
    AND S.record_status != 'D'
    $sample_category_clause
    -- ORDER BY PISB.n_observations, $args{peptide_field} ASC
    $order_by 
  ~;
  my @samples = $sbeams->selectSeveralColumns($sql);

  my $sample_js;
  my %samples;
  my %peptides;
  my $cntr = 0;
  for my $row ( @samples ) { 
    $cntr++;
    my $key = $row->[1] . '::::' . $row->[0];
    $row->[3] .= '*' if $row->[4] < 2;
    $peptides{$row->[3]}++;
    $samples{$key} ||= {};
    $samples{$key}->{$row->[3]} += $row->[2];
  }
  my $array_def = qq~
    <script type="text/javascript">
    google.setOnLoadCallback(drawHeatMap);
    function drawHeatMap() {
    var data = new google.visualization.DataTable();
    data.addColumn('string', 'Sample Name');
    ~;
  my @peps = ( $args{force_order} ) ? @{$args{force_order}} : sort( keys( %peptides ) );
#  for my $pa ( sort( keys( %peptides ) ) ) {
  for my $pa ( @peps ) {
    $array_def .= "    data.addColumn('number', '$pa');\n";
  }
  $array_def .= "DEFINE_ROWS_HERE\n";

  my $row = 0;
  my $max = 0;
  my $min = 50;
  $cntr = 0;
  for my $sa ( sort { "\L$a" cmp "\L$b" } ( keys( %samples ) ) ) {
    my $col = 0;
    my ( $name, $id ) = split "::::", $sa;
    $array_def .= "    data.setValue( $row, $col, '$name' );\n";
    $col++;

    for my $pa ( @peps ) {
      if ( $samples{$sa}->{$pa} ) {
        my $pep_cnt = log(1 + $samples{$sa}->{$pa})/log(10);
	$max = ( $pep_cnt > $max ) ? $pep_cnt : $max;
	$min = ( $pep_cnt < $min ) ? $pep_cnt : $min;
        $array_def .= "    data.setValue( $row, $col, $pep_cnt );\n";
	$cntr++;
      }
      $col++;
    }
    $row++;
  }
  $array_def =~ s/DEFINE_ROWS_HERE/data.addRows($row);/;
  my $num_colors = 256;
  $array_def .= qq~
    heatmap = new org.systemsbiology.visualization.BioHeatMap(document.getElementById('heatmapContainer'));
    heatmap.draw(data, {numberOfColors:$num_colors,passThroughBlack:false,startColor:{r:255,g:255,b:255,a:1},endColor:{r:100,g:100,b:100,a:1},emptyDataColor:{r:256,g:256,b:256,a:1}});
		}
  </script>
  ~;

  $args{header_text} = ( $args{header_text} ) ? "<TR $trinfo><TD ALIGN=CENTER CLASS=section_description>$args{header_text}</TD></TR>" : '';
  $args{second_header} = ( $args{second_header} ) ? "<TR $trinfo><TD ALIGN=CENTER CLASS=info_text>$args{second_header}</TD></TR>" : '';
  my $content = qq~
  <script type="text/javascript" src="https://ajax.googleapis.com/ajax/libs/prototype/1.6.1.0/prototype.js"></script>
  <script type="text/javascript" src="https://www.google.com/jsapi"></script>
  <script type="text/javascript">
    google.load("visualization", "1", {});
  </script>
  <script type="text/javascript" src="$HTML_BASE_DIR/usr/javascript/main/js/load.js"></script>
  <script type="text/javascript">
    systemsbiology.load("visualization", "1.0", {packages:["bioheatmap"]});
  </script>
  $array_def

  $args{header_text}
  $args{second_header}
	<TR $trinfo><TD> <DIV ID="heatmapContainer"></DIV>  </TD></TR>
	~;

  return ( wantarray() ) ? ($header, $content) : $header . "\n" . $content;

} # end getSampleMapDisplay

sub getBuildSelector {
  my $self = shift;
  my %args = @_;
  my $build_id = $args{atlas_build_id};
  my $sbeams = $self->getSBEAMS();

  my $accessible_builds = join( ',', $self->getAccessibleBuilds() );
  my $accessible_projects = join( ',', $sbeams->getAccessibleProjects() );

  # Get a hash of available atlas builds
  my $sql = qq~
  SELECT atlas_build_id, atlas_build_name
  FROM $TBAT_ATLAS_BUILD
  WHERE project_id IN ( $accessible_projects )
  AND record_status!='D'
  ORDER BY atlas_build_name
  ~;
  my @ordered_ids;
  my %id2build;
  my $sth = $sbeams->get_statement_handle( $sql );
  while ( my @row = $sth->fetchrow_array() ) {
    push @ordered_ids, $row[0];
    $id2build{$row[0]} = $row[1];
    $build_id ||= $row[0];
  }
  my $build_selector =  $q->popup_menu( -name => "atlas_build_id",
                                        -values => [ @ordered_ids ],
                                        -labels => \%id2build,
                                        -default => $build_id,
                                        -onChange => 'switchAtlasBuild()' );

  my $style = '';
  if ($args{inline}) {
    $style .= 'display:inline;';
  }
  if ($style) {
    $style = "style='$style'";
  }

  my $selector_widget = qq~
    <form $style name="build_form" id="build_form">
    $build_selector
    </form>
    <script LANGUAGE="Javascript">
      function switchAtlasBuild() {
        document.build_form.submit();
      }
    </script>
  ~;
  return $selector_widget;
}

sub getSampleMapDisplayMod {
  my $self = shift;
  my $sbeams = $self->getSBEAMS();
  my %args = ( peptide_field => 'peptide_accession', @_ );

  my @all_peps = ( @{$args{snp_support}}, @{$args{snp_original}} );
  my @speps = @{$args{snp_support}};
  my $in = join( ", ", @all_peps );
  return unless $in;

  my $header = '';
  if ( $args{link} ) {
    $header .= $self->encodeSectionHeader( text => 'Sample peptide map:',
                                          link => $args{link},
                                          anchor => 'samples'
                                          );
  } else {
    $header .= $self->encodeSectionHeader( text => 'Sample peptide map:',
                                          anchor => 'samples'
                                         );
  }
  $header = '' if $args{no_header};

  my $html = '';
  my $trinfo = $args{tr_info} || '';


  my $sql = qq~     
  	SELECT DISTINCT SB.atlas_search_batch_id, sample_tag, 
		PISB.n_observations, peptide_sequence, PI.peptide_instance_id
		FROM $TBAT_ATLAS_SEARCH_BATCH SB 
	  JOIN $TBAT_SAMPLE S ON s.sample_id = SB.sample_id
	  JOIN $TBAT_PEPTIDE_INSTANCE_SEARCH_BATCH PISB ON PISB.atlas_search_batch_id = SB.atlas_search_batch_id
	  JOIN $TBAT_PEPTIDE_INSTANCE PI ON PI.peptide_instance_id = PISB.peptide_instance_id
	  JOIN $TBAT_PEPTIDE P ON P.peptide_id = PI.peptide_id
    WHERE PI.peptide_instance_id IN ( $in )
    AND S.record_status != 'D'
    ORDER BY sample_tag ASC
  ~;

  my @samples = $sbeams->selectSeveralColumns($sql);

  my $sample_js;
	my %samples;
	my %peptides;
	my $cntr = 0;
  my %id2seq;
  for my $row ( @samples ) { 
		$cntr++;
		my $key = $row->[1] . '::::' . $row->[0];
		$peptides{$row->[3]}++;
		$samples{$key} ||= {};
		$samples{$key}->{$row->[3]} = $row->[2];
    $id2seq{$row->[4]} = $row->[3];
	}
	my $array_def = qq~
	<script type="text/javascript">
    google.setOnLoadCallback(drawHeatMap);
    function drawHeatMap() {
    var data = new google.visualization.DataTable();
    data.addColumn('string', 'Sample Name');
	~;
#	my @peps = ( $args{force_order} ) ? @{$args{force_order}} : sort( keys( %peptides ) );
  my %speps;
  for my $pep ( @speps ) {
    $speps{$pep}++;
  }

#  for my $pa ( sort( keys( %peptides ) ) ) {
  for my $pi ( @all_peps ) {
    my $pa = $id2seq{$pi};
    $array_def .= "    data.addColumn('number', '$id2seq{$pi}');\n";
	}
	$array_def .= "DEFINE_ROWS_HERE\n";

	my $row = 0;
	my $max = 0;
	my $min = 50;
	$cntr = 0;
	for my $sa ( sort( keys( %samples ) ) ) {
	  my $col = 0;
		my ( $name, $id ) = split "::::", $sa;
    $array_def .= "    data.setValue( $row, $col, '$name' );\n";
	  $col++;

    for my $pi ( @all_peps ) {
      my $pa = $id2seq{$pi};
			if ( $samples{$sa}->{$pa} ) {
        my $pep_cnt = log(1 + $samples{$sa}->{$pa})/log(10);
			  $max = ( $pep_cnt > $max ) ? $pep_cnt : $max;
    		$min = ( $pep_cnt < $min ) ? $pep_cnt : $min;
        $array_def .= "    data.setValue( $row, $col, $pep_cnt );\n";
		    $cntr++;
			}
		  $col++;
		}
		$row++;
	}
	$array_def =~ s/DEFINE_ROWS_HERE/data.addRows($row);/;
	my $num_colors = 256;
	$array_def .= qq~
	heatmap = new org.systemsbiology.visualization.BioHeatMap(document.getElementById('heatmapContainer'));
	heatmap.draw(data, {numberOfColors:$num_colors,passThroughBlack:false,startColor:{r:255,g:255,b:255,a:1},endColor:{r:100,g:100,b:100,a:1},emptyDataColor:{r:256,g:256,b:256,a:1}});
		}
  </script>
  ~;

	$args{header_text} = ( $args{header_text} ) ? "<TR $trinfo><TD ALIGN=CENTER CLASS=section_description>$args{header_text}</TD></TR>" : '';
	$args{second_header} = ( $args{second_header} ) ? "<TR $trinfo><TD ALIGN=CENTER CLASS=info_text>$args{second_header}</TD></TR>" : '';
	my $content = qq~
  <script type="text/javascript" src="https://ajax.googleapis.com/ajax/libs/prototype/1.6.1.0/prototype.js"></script>
  <script type="text/javascript" src="https://www.google.com/jsapi"></script>
  <script type="text/javascript">
    google.load("visualization", "1", {});
  </script>
  <script type="text/javascript" src="$HTML_BASE_DIR/usr/javascript/main/js/load.js"></script>
  <script type="text/javascript">
    systemsbiology.load("visualization", "1.0", {packages:["bioheatmap"]});
  </script>

	$array_def

  $args{header_text}
  $args{second_header}
	<TR $trinfo><TD> <DIV ID="heatmapContainer"></DIV>  </TD></TR>
	~;

  return ( wantarray() ) ? ($header, $content) : $header . "\n" . $content;


} # end getSampleMapDisplay



### END

###############################################################################
# display enhance sample info
###############################################################################
sub getPeptideSampleDisplay {
  my $self = shift;
  my $sbeams = $self->getSBEAMS();
  my %args = @_;
  my $SUB_NAME = 'getDetailedSampleDisplay';
  my $rows_to_show = $args{rows_to_show} || 25;
  my $mia = '';
  for my $arg ( qw ( sample_ids build_clause peptide_clause ) ) {
    next if defined $args{$arg};
    my $sep = ( $mia ) ? ',' : '';
    $mia .= $mia . $sep . $arg;
  }
  $args{tr_info} ||= '';
  if ( $mia ) {
    $log->error( "Missing required argument(s} $mia" );
    return;
  }
  my $in = join( ", ", @{$args{sample_ids}} );
  return unless $in;

  my $sql = qq~
  	SELECT S.sample_id, sample_tag, PISB.n_observations,
           instrument_name, CASE WHEN ENZ.name IS NULL THEN 'Trypsin' ELSE ENZ.name END AS Enzyme,
           PUB.publication_name, PUB.abstract , PUB.uri, S.repository_identifiers
		FROM $TBAT_ATLAS_SEARCH_BATCH SB 
	  JOIN $TBAT_SAMPLE S ON s.sample_id = SB.sample_id
	  JOIN $TBAT_PEPTIDE_INSTANCE_SEARCH_BATCH PISB ON PISB.atlas_search_batch_id = SB.atlas_search_batch_id
	  JOIN $TBAT_PEPTIDE_INSTANCE PI ON PI.peptide_instance_id = PISB.peptide_instance_id
	  JOIN $TBAT_PEPTIDE P ON P.peptide_id = PI.peptide_id
	  LEFT JOIN $TBPR_INSTRUMENT I ON S.instrument_model_id = I.instrument_id 
	  LEFT JOIN $TBAT_PROTEASES ENZ ON ENZ.id = S.protease_id
    LEFT JOIN $TBAT_SAMPLE_PUBLICATION SP ON (SP.sample_id = S.sample_id and SP.record_status != 'D') 
    LEFT JOIN $TBAT_PUBLICATION PUB ON PUB.publication_id = SP.publication_id  
    WHERE S.sample_id IN ( $in )
    $args{build_clause}
    $args{peptide_clause}
    AND S.record_status != 'D'
    ORDER BY S.sample_id
  ~;
  my @rows = $sbeams->selectSeveralColumns($sql);
  my $table = $self -> getSampleTableDisplay (data => \@rows,
                  rows_to_show => $rows_to_show,
                  type => 'Peptide');
  return $table;
}


sub getSampleTableDisplay{
  my $self = shift;
  my %args = @_;
  my $data = $args{data};
  my $type = $args{type};
  my $rows_to_show = $args{rows_to_show} || 25;
  my $sample_category = $args{sample_category} || 0;
  my @samples = ();
  my $pre_id = '';
  my $anno =0;
  my @annotation_urls;
  foreach my $row (@$data){
     my ($id, $title, $nobs,$ins,$enzyme,$pub_name, $abstract, $link, $rp_id,$sc) = @$row;
     my $annotation_url;
     $id = "<a href='". "$CGI_BASE_DIR/PeptideAtlas/ManageTable.cgi?TABLE_NAME=AT_SAMPLE&sample_id=" . $id . "' target='_blank'>$id</a>";
     ($rp_id, $annotation_url) = $self->get_dataset_url ($rp_id);
     $anno = 1 if ($annotation_url ne '');
     push @annotation_urls, $annotation_url;
     if ($pre_id eq $id){next;$pre_id = $id;}  ## keep one publication only
     if ($abstract){
       $pub_name = $self->make_pa_tooltip( tip_text => "Abstract: $abstract", link_text => "<a href='$link'>$pub_name</a>" );
     }
     if ($type eq 'Peptide'){
       push @samples, [$id, $rp_id, $title, $nobs,$ins,$enzyme,$pub_name];
     }else{
       if ($sample_category){
				 push @samples, [$id, $rp_id, $title,$ins,$enzyme,$pub_name,$sc];
       }else{
         push @samples, [$id, $rp_id, $title,$ins,$enzyme,$pub_name]; 
       }
     }
     $pre_id = $id;
  }
  my @align=();
  my @headings;
  my @cols = ();
  if ($type eq 'Peptide'){
    @cols = ('Experiment ID','Dataset','Experiment Tag','NObs','Instrument','Enzyme','Publication');
    @align = qw(center left left center left left);
  }else{
    if ($sample_category){
      @cols = ('Experiment ID','Dataset','Experiment Tag', 'Instrument','Enzyme','Publication', 'Source');
    }else{
      @cols = ('Experiment ID','Dataset','Experiment Tag', 'Instrument','Enzyme','Publication');
    }
    @align = qw(center left left center left left left);
  }
  if ($anno){
    for (my $i=0; $i<=$#samples; $i++){
      splice @{$samples[$i]}, 2, 0, $annotation_urls[$i];
    }
    splice @cols, 2, 0, 'Experiment Annotation';
    splice @align, 2, 0, 'center';

  }  

  foreach my $col (@cols){
    push @headings, $col, '';
  }
  my $headings = $self->make_sort_headings( headings => \@headings,  default => 'Experiment ID');
  unshift @samples, ($headings);
  my $table = $self->encodeSectionTable( header => 1, 
					 header_sticky => 1,
                                         tr_info => $args{tr_info},
                                         align  => [@align],
					 bkg_interval => 3,
                                         nowrap => [qw(4 6)],
                                         rows_to_show => $rows_to_show,
                                         max_rows => $args{max_rows},
                                         sortable => 1,
                                         rows => \@samples );
  return $table;
}


sub getPTMTableDisplay {
  my $self = shift;
  my %args = @_;
  my $SUB_NAME = 'getPTMTableDisplay';
  my $cols = $args{cols};
  my $data = $args{data};
  my $ptm_type = $args{ptm_type};
  my $atlas_build_id = $args{atlas_build_id}; 
  my $biosequence_name = $args{biosequence_name};
  my @rows;
  shift @$cols;
  pop @$cols;
  push @rows, [@$cols];
  my $n_cols = scalar @$cols;
  foreach my $prot (keys %$data){
    foreach my $pos (sort {$a <=> $b} keys %{$data->{$prot}}){
      foreach my $vals (@{$data->{$prot}{$pos}}){
				my @row = ();
        my $site = $pos+1;
        my $residue=$vals->[0]; 
				push @row,$site; 
        my $url="$CGI_BASE_DIR/PeptideAtlas/GetPTMSpectra?".
                  "atlas_build_id=$atlas_build_id&site=$site&biosequence_name=$prot".
                  "&ptm_type=$ptm_type&residue=$residue&apply_action=QUERY";

        #my @columns = ('Residue','FLR Category', 'nObs', 'One_mod', 'Two_mods', 'Over_two_mods',
        #         'nP<.01', 'nP<.05', 'nP<.20', 'nP.2-.8', 'nP>.80', 'nP>.95', 'nP>.99',
        #                  'no-choice','enriched-with-mod','enriched-but-non-mod','non-enriched',
        #                           'InNextProt','InUniprot','peptide');
				if ($vals->[2] > 0 || $vals->[17] eq 'yes' || $vals->[18] eq 'yes'){
					my $start_in_biosequence = $pos + 1;
					my $link = "$CGI_BASE_DIR/PeptideAtlas/GetPeptide?".
										 "atlas_build_id=$atlas_build_id&searchWithinThis=Peptide+Sequence&searchForThis=".
										 "$vals->[19]&apply_action=QUERY"; 
					$vals->[0] = $self->make_pa_tooltip( tip_text => "Get most observed peptide sequence covering this site",
																								link_text => "<a href='$link' target='_blank'>$vals->[0]</a>" );
				}
				foreach my $i (0..$n_cols-2){ 
					if ( $vals->[$i] eq '0' ||  $vals->[$i] eq 'no'){
						$vals->[$i] = '-';
					}
          ## make links for nP>.95 and nP>.99 ids 
          my $link = $url; 
          if ($vals->[10] > 0){
            $link = $url ."&min=0.80&max=0.95";
            $vals->[10] = $self->make_pa_tooltip( tip_text => "Get observed spectra covering this site with ptm scores within 0.80-0.95",
                                                 link_text => "<a href='$link' target='_blank'>$vals->[10]</a>" );
          }
          if ($vals->[11] > 0){
            $link =$url ."&min=0.95&max=0.99";
            $vals->[11] = $self->make_pa_tooltip( tip_text => "Get observed spectra covering this site with ptm scores within 0.95-0.99",
                                                 link_text => "<a href='$link' target='_blank'>$vals->[11]</a>" );
          }
          if ($vals->[12] > 0){
            $link = $url ."&min=0.99";
            $vals->[12] = $self->make_pa_tooltip( tip_text => "Get observed spectra covering this site with ptm scores within 0.99-1.0",
                                                 link_text => "<a href='$link' target='_blank'>$vals->[12]</a>" );
          }
          if ($vals->[13] > 0){
            $link = $url ."&nochoice=1";
            $vals->[13] = $self->make_pa_tooltip( tip_text => "Get observed spectra covering this site that have no choice in the localization of the PTM ", link_text => "<a href='$link' target='_blank'>$vals->[13]</a>" );
          }
          push@row, $vals->[$i]; 
			  }	
        push @rows , [@row];
      }
    }
  }
  my @align = ();
  foreach my $i(0..15){
    push @align, 'center';
  }
  my $tr_info = '';
  my $table = '';
  if (@rows > 15){
    $table = $self->encodeSectionTable( header => 1,
					header_sticky => 1,
					y_scroll => 'style="overflow-y:scroll; height:300px; display:block;" ', 
					unified_widgets => 1,
					set_download => 1,
					align  => [@align],
					bkg_interval => 3,
					file_prefix => 'ptm_',
					rows => \@rows );
  }else{
    $table = $self->encodeSectionTable( header => 1,
					unified_widgets => 1,
					set_download => 1,
					align  => [@align],
					bkg_interval => 3,
					file_prefix => 'ptm_',
					rows => \@rows );
    
  }
  return $table;
}


###############################################################################
# get_individual_spectra_display 
###############################################################################
sub get_individual_spectra_display {
  my $self = shift;
  my $sbeams = $self->getSBEAMS();
  my %args = @_;
  my $SUB_NAME = 'get_individual_spectra_display';
  my $bg_color = $args{bg_color} || '';
  my $column_titles_ref = $args{column_titles_ref};
  my $colnameidx_ref = $args{colnameidx_ref};
  my $hidden_cols_ref = $args{hidden_cols_ref};
  my $table_id = $args{table_id} || 'individual_spectra'; 

  my $sortable = 1;
  if ( $args{sortable} ne ''){
    $sortable = $args{sortable};
  }


  unless( $args{resultset_ref} ) {
    $log->error( "No resultset_ref passed to display" );
    return;
  }
  my $resultset_ref = $args{resultset_ref};
  my @data = ();
  my @headings =();

  my $loop =1;

  foreach my $row (@{$resultset_ref->{data_ref}}){
    my @tmp =();
    foreach my $col (sort {$colnameidx_ref->{$a} <=> $colnameidx_ref->{$b}} keys %$colnameidx_ref){
       if (! $hidden_cols_ref->{$col}){
         push @tmp, $row->[$colnameidx_ref->{$col}];
         push @headings, $column_titles_ref->[$colnameidx_ref->{$col}] if ($loop==1);
       }
    }
    $loop++;
    push @data, \@tmp;
  }
  my $headings = $self->get_column_defs( labels => \@headings, plain_hash => 1 );
  unshift @data, ($self->make_sort_headings(headings => $headings));

  
    #push @{$resultset_ref->{column_list_ref}}, 'num_prot_mappings';
    #push @{$resultset_ref->{types_list_ref}}, 'int';

  my $align = [qw(left center center center left center center center center center center center center center left center center)];


  my $html = '';
  if (@data > 15){
     $html =  $self->encodeSectionTable( header => 1,
          header_sticky => 1,
          y_scroll => 'style="overflow-y:scroll; height:300px; display:block;" ',
					unified_widgets => 1,
					set_download => 1,
					colspan => 3,
					align  => $align,
					rows => \@data,
					max_rows => 500,
					nowrap => [1],
					bkg_interval => 3,
					rs_headings => \@$column_titles_ref, 
					bg_color => '#EAEAEA',
					sortable => 1,
					close_table => 1,
					table_id => $table_id,
					truncate_msg => '  <b>Use link in Modified Peptides section to filter results</b>.',
      );
    }else{
     $html =  $self->encodeSectionTable( header => 1,
					unified_widgets => 1,
					set_download => 1,
					colspan => 3,
					align  => $align,
					rows => \@data,
					nowrap => [1],
					bkg_interval => 3,
					rs_headings => \@$column_titles_ref, 
					bg_color => '#EAEAEA',
					sortable => 1,
					close_table => 1,
					table_id => $table_id,
      );
    }
    return "<TABLE>$html\n";


} # end 

###############################################################################
# displaySamples
###############################################################################
sub getProteinSampleDisplay {
  my $self = shift;
  my $sbeams = $self->getSBEAMS();
  my %args = @_;
  my $SUB_NAME = 'getProteinSampleDisplay';
  my $bg_color = $args{bg_color} || '';
  my $sortable = 1;
  if ( $args{sortable} ne ''){
    $sortable = $args{sortable};
  }
  my $rows_to_show = $args{rows_to_show} || 25;
  my $sample_category = $args{sample_category} || 0;
  unless( $args{sample_ids} ) {
    $log->error( "No samples passed to display samples" );
    return;
  }

  my $in = join( ", ", @{$args{sample_ids}} );
  return unless $in;
	my $sql = qq~
	SELECT S.SAMPLE_ID,
				 S.SAMPLE_Tag,
				 '' ,
				 INSTRUMENT_NAME,
				 CASE WHEN ENZ.name IS NULL THEN 'Trypsin' ELSE ENZ.name END AS Enzyme,
				 PUB.PUBLICATION_NAME,
				 PUB.ABSTRACT ,
				 PUB.URI,
				 S.REPOSITORY_IDENTIFIERS,
				 SC.name 
	FROM $TBAT_SAMPLE S
	LEFT JOIN $TBPR_INSTRUMENT I ON S.instrument_model_id = I.instrument_id
	LEFT JOIN $TBAT_PROTEASES ENZ ON ENZ.id = S.protease_id
	LEFT JOIN $TBAT_SAMPLE_PUBLICATION SP ON SP.SAMPLE_ID = S.SAMPLE_ID
	LEFT JOIN $TBAT_PUBLICATION PUB ON (PUB.PUBLICATION_ID = SP.PUBLICATION_ID
																			AND SP.record_status != 'D')
	JOIN $TBAT_SAMPLE_CATEGORY SC ON (SC.id = S.sample_category_id)
	WHERE S.SAMPLE_ID IN ( $in )
	AND S.RECORD_STATUS != 'D'
	ORDER BY sample_ID
	~;

  my @rows = $sbeams->selectSeveralColumns($sql);
  my $html = ''; 
  if ($sample_category){
    ## get pie chart 
    my %sc = ();
    foreach my $row(@rows){
      $sc{$row->[9]}++;
    }
    my @sc = keys %sc;
    my @vals = values %sc;
		$html = qq~
			<div  id ='sample_cat_pie' style='display:block'>
		 ~;
		$html .= $self->plotly_pie(data => \@vals, names=>\@sc , divName=>"sample_cat_pie"); 
  }
  $html .= $self -> getSampleTableDisplay(data => \@rows,
                               rows_to_show => $rows_to_show, 
                               type => 'Protein',
                               sample_category=>$sample_category,
                               );
  return $html;
} # end getProteinSampleDisplay 

sub add_tabletoggle_js {
  my $self = shift;
  return '' if $self->{_added_ttoggle_js};
  $self->{_added_ttoggle_js}++;
  return <<"  END";
  <STYLE TYPE="text/css" media="screen">
    tr.visible { display: block-row; }
    tr.hidden { display: none; }
    td.visible { display: table-cell; }
    td.hidden { display: none; }
  </STYLE>
  <SCRIPT TYPE="text/javascript">
    function toggle_em(prefix) {
      
      // Grab page elements by their IDs
      var togglekey = prefix + '_toggle';
      var textkey = prefix   + '_text';

      var rows = document.getElementsByName(togglekey);
      var show = document.getElementById(textkey);

      var ttext = show.innerHTML;
      if ( ttext == '[Show more rows]' ) {
        show.innerHTML = "[Show fewer rows]";
      } else {
        show.innerHTML = "[Show more rows]";
      }
      
      for (var i=0; i < rows.length; i++) {
        if ( rows[i].className == 'hidden' ) {
           rows[i].className = 'visible';
        } else {
           rows[i].className = 'hidden';
        }
      }
    }
    </SCRIPT>

  END
  
}

#+
# Routine will set up PA tooltip colors and javascript.  
#
# @return   JS/CSS to be printed in the page at an appropriate juncture
#-
sub init_pa_tooltip {
  my $self = shift;
  
  # return nothing if tooltip has already been initialized
  return '' if $self->{_tooltip_init};
  
  $self->{_tooltip_init} = 1;
  return <<"  END";
  <SCRIPT DEFER="defer" src="$HTML_BASE_DIR/usr/javascript/TipWidget.js" TYPE="text/javascript"></SCRIPT>
  <STYLE>
  div#tooltipID { background-color:#F0F0F0;
                  border:2px 
                  solid #FF8C00; 
                  padding:4px; 
                  line-height:1.5; 
                  width:auto; 
                  font-family: Helvetica, Arial, sans-serif;
                  font-size:12px; 
                  font-weight: normal; 
                  position:absolute; 
                  visibility:hidden; left:0; top:0;
                }
  </STYLE>
  END
}


#+
# Routine to generate a 'tooltip' mouseover widget with 
#-
sub make_pa_tooltip {
  my $self = shift;
  my %args = @_;
  for my $arg ( qw( tip_text link_text ) ) {
    if ( !$args{$arg} ) {
      $log->error( "Missing required argument $arg" );
      return "";
    }
  }

  # Issue warning if tooltip js wasn't printed (well, fetched anyway).
  $self->init_pa_tooltip();
  $log->warn( "Failed to itialize tooltip!" ) if !$self->{_tooltip_init};

  my $class = $args{class} || 'pseudo_link';

  # Sanitize tooltip text
  $args{tip_text} =~ s/\'/\\\'/g;
  $args{tip_text} =~ s/\"/\\\'/g;
  $args{tip_text} =~ s/\r\n/ /g;
  $args{tip_text} =~ s/\n/ /g;
  
  return "<SPAN CLASS=$class onMouseover=\"showTooltip(event, '$args{tip_text}')\" onMouseout=\"hideTooltip()\">$args{link_text}</SPAN>";
}

sub formatMassMods { 
  my $self = shift;
  my $sequence = shift || return undef;
  $sequence =~ s/([\(\[])/<SPAN CLASS="aa_mod">$1/gm;
  $sequence =~ s/([\]\)])/$1<\/SPAN>/gm;
  return $sequence;
}


sub get_table_help_section {
  my $self = shift;
  my %args = @_;
  $args{showtext} ||= 'show column descriptions';
  $args{hidetext} ||= 'hide column descriptions';
  $args{heading} ||= '';
  $args{description} ||= '';
  $args{footnote} ||= '';


  my $ecnt = 0;
  my $index = "<TABLE class=info_box>\n";
  for my $entry ( @{$args{entries}} ) {
    $ecnt++;
    $index .= $self->encodeSectionItem( %$entry, nowrap_description => 1 );
  }
  $index .= "</TABLE>\n";

  my $content =<<"  END";
  <BR>
  <span class=section_heading>$args{heading}</span> 
  <span class=description>$args{description}</span>
  $index
  <span class=description>$args{footnote}</span>
  END

  $sbeams = $self->getSBEAMS();
  my $section_toggle = $sbeams->make_toggle_section( content => $content,
                                                     sticky   => 0,
                                                     visible  => 0,
                                                     imglink  => 1,
                                                     showimg  => "/info_small.gif",
                                                     hideimg  => "/info_small.gif",
                                                     textlink => 1,
                                                     name     => $args{name},
                                                     showtext => $args{showtext},
                                                     hidetext => $args{hidetext},
                                          );

  return "$section_toggle<BR>";
}


sub get_table_help {
  my $self = shift;
  my %args = @_;
  my $column_titles_ref = $args{column_titles_ref};
  my $colnameidx_ref = $args{colnameidx_ref};
  my $hidden_cols_ref = $args{hidden_cols_ref};
  $args{heading} ||= '';
  $args{description} ||= '';
  $args{footnote} ||= '';
  $args{showtext} ||= 'show column descriptions';
  $args{hidetext} ||= 'hide column descriptions';


  my @headings = ();
  if ($colnameidx_ref && $hidden_cols_ref){ 
    foreach my $col (sort {$colnameidx_ref->{$a} <=> $colnameidx_ref->{$b}} keys %$colnameidx_ref){
      if (not defined $hidden_cols_ref->{$col}){
	push @headings, $column_titles_ref->[$colnameidx_ref->{$col}];
      }
    }
  }else{
   @headings = @$column_titles_ref;
  }

  my $headings = $self->get_column_defs( labels => \@headings );

  my $index = "<TABLE class=info_box>\n";
  for my $entry ( @$headings ) {
    $index .= $self->encodeSectionItem( %$entry, nowrap_description => 1 );
  }
  $index .= "</TABLE>\n";

  my $content =<<"  END";
  <BR>
  <span class=section_heading>$args{heading}</span> 
  <span class=description>$args{description}</span>
  $index
  <span class=description>$args{footnote}</span>
  END

  $sbeams = $self->getSBEAMS();
  my $section_toggle = $sbeams->make_toggle_section( content => $content,
                                                     sticky   => 0,
                                                     visible  => 0,
                                                     imglink  => 1,
                                                     showimg  => "/info_small.gif",
                                                     hideimg  => "/info_small.gif",
                                                     textlink => 1,
                                                     name     => $args{name},
                                                     showtext => $args{showtext},
                                                     hidetext => $args{hidetext},
                                          );

  return "$section_toggle<BR>";
}

###############################################################################
########### vocabHTML
######################################################################################

sub vocabHTML {

###### This subroutine displays the Vocab section on the HTML Page
 
my $self=shift;
my $sbeamsMOD=$self;

my ($table3);

my $pi_desc='Isoelectric point of the peptide';
my $ssrcalc_desc='Sequence Specific Retention Factor provides a hydrophobicity measure for each peptide using the algorithm of Krohkin et al. Version 3.0 <A HREF=http://hs2.proteome.ca/SSRCalc/SSRCalc.html target=_blank>more</A>';
my $org_desc='Organism in which this peptide is observed';
my $n_obs_desc='Number of MS/MS spectra that are identified to this peptide in each build. The hyperlink will take you to the peptide page that will display all the relevant information about this peptide contained in the listed PeptideAtlas build';
my $build_names_desc='Public build in which this peptide is found. The hyperlink will take you to a build specific page summarizing all the relevant information about the build.';

$table3 ='<BR><BR><BR><BR><BR><BR><BR><BR><BR><BR><TABLE WIDTH=600>';

$table3 .= $sbeamsMOD->encodeSectionHeader( text => 'Vocabulary',
                                            width => 900
                                            );

#Substituting the path generated from encodeSectionHeader to match with peptideatlas.org for displaying orange header
$table3 =~ s/\/devAP\/sbeams//gm;

$table3 .= $sbeamsMOD->encodeSectionItem( key   => 'pI',
                                          value => $pi_desc,
                                          key_width => '20%'
                                         );
$table3.= $sbeamsMOD->encodeSectionItem( key   => 'SSRCalc',
                                          value =>$ssrcalc_desc
                                         );
$table3.= $sbeamsMOD->encodeSectionItem( key   => 'Organism Name',
                                          value =>$org_desc
                                         );

$table3.= $sbeamsMOD->encodeSectionItem( key   => 'Number of observations ',
                                          value =>$n_obs_desc
                                         );

$table3.= $sbeamsMOD->encodeSectionItem( key   => 'Build Names in which Peptide is found',
                                          value =>$build_names_desc
                                         );
$table3 .= '</TABLE>';

return ($table3);


}





sub get_atlas_checklist {
  my $self = shift;
  my %args = @_;

  # check input
	for my $arg ( qw( build_info js_call ) ) {
		return "Missing required arguement $arg" unless defined $args{$arg};
	}

  #### Get a list of accessible builds, use to validate passed list.
  my @accessible = $self->getAccessibleBuilds();
	my %build_check;
	for my $build ( @accessible ) {
		$build_check{$build}++;
	}

  # Convenience
	my %build_info = %{$args{build_info}};
  my @builds =  sort { $build_info{$a}->{org} cmp $build_info{$b}->{org} ||
		               $build_info{$b}->{is_curr} <=>  $build_info{$a}->{is_curr} ||
		            $build_info{$b}->{is_default} <=>  $build_info{$a}->{is_default} ||
						 		      $build_info{$a}->{name} cmp $build_info{$b}->{name} } keys %build_info;

  my $table = SBEAMS::Connection::DataTable->new(BORDER => 0);
  $table->addRow( [ 'Add', 'Build Id', 'Display Name', 'Full Name', 'Organism', 'is_def' ] );
  $table->setRowAttr(  ROWS => [1], BGCOLOR => '#bbbbbb', ALIGN=>'CENTER' );
  $table->setHeaderAttr( BOLD => 1 );

  for my $build ( @builds ) {
		my %build = %{$build_info{$build}};
    my $checked = ( $build{visible} ) ? 'checked' : '';

    my $chkbox =<<"    END";
    <INPUT $checked TYPE="checkbox" NAME="build_id" VALUE="$build" onchange="$args{js_call}($build);">
    END
    $table->addRow( [ $chkbox, $build, @build{qw(display_name name org is_default)} ] );
    $table->setRowAttr( ROWS => [$table->getRowNum()], BGCOLOR => $build{bgcolor} );

	}
	return "$table";
}

sub get_atlas_select {

  my $self = shift;
  my %args = @_;
  my $select_name = $args{select_name} || 'build_id';

  #### Get a list of accessible builds, use to validate passed list.
  my @accessible = $self->getAccessibleBuilds();
  my $build_string = join( ",", @accessible );
  return '' unless $build_string;

  my $sql = qq~
  SELECT atlas_build_id, build_tag 
  FROM $TBAT_ATLAS_BUILD
  WHERE atlas_build_id IN ( $build_string )
  ~;

  $log->info( $sql );

  my $select = "<INPUT TYPE=SELECT NAME=$select_name>\n";
  my $sth = $sbeams->get_statement_handle->( $sql );
  while ( my @row = $sth->fetchrow_array() ) {
    $select .= "<OPTION VALUE=$row[0]>$row[1]\n";
  }
  $select .= "</INPUT>\n";
  return $select;
}

sub get_proteome_coverage_new {
  my $self = shift;
  my $build_id = shift; 
  my $patterns = shift;
  my $data_only = shift;
  my @patterns = @$patterns;

  return [] if (! $build_id || ! @patterns);
  my $sbeams = $self->getSBEAMS();
  my $sql = '';
  my $obs_sql = '';
  my @names = ();
  my @name_desc =();
  my %biosqeuences = ();
  my %entry_cnt = ();
  my %obs=();
  my %ptm_summary =();
  my $sql = qq~
		  SELECT B.BIOSEQUENCE_ID,B.BIOSEQUENCE_NAME, B.BIOSEQUENCE_Desc, B.BIOSEQUENCE_seq
			FROM $TBAT_ATLAS_BUILD AB
			JOIN $TBAT_BIOSEQUENCE B ON B.BIOSEQUENCE_SET_ID = AB.BIOSEQUENCE_SET_ID
			WHERE atlas_build_id = $build_id
  ~;
  my @rows = $sbeams->selectSeveralColumns($sql);
  foreach my $row(@rows){
    my ($id, $name, $desc,$seq) = @$row;
    $biosqeuences{$id}{accession} = $name;
    $biosqeuences{$id}{desc} = $desc;
    $biosqeuences{$id}{seq} = $seq;
  }

  my $obs_sql = qq~
		SELECT DISTINCT B2.BIOSEQUENCE_ID, B2.BIOSEQUENCE_ID
		FROM $TBAT_PEPTIDE_INSTANCE PI
		JOIN $TBAT_PEPTIDE_MAPPING PM ON (PI.PEPTIDE_INSTANCE_ID = PM.PEPTIDE_INSTANCE_ID)
		JOIN $TBAT_BIOSEQUENCE B2 ON (B2.biosequence_id = PM.matched_biosequence_id)
		WHERE PI.atlas_build_id = $build_id
		--AND B2.BIOSEQUENCE_NAME LIKE 'CONTAM%' 
   ~;
  my %biosqeuence_id_obs = $sbeams->selectTwoColumnHash($obs_sql);
 
  foreach my $line (@patterns){
    my ($org_id, $name, $type, $pat_str, $pat_desc)  = split(/\t/, $line);
    my @pats = split(/;/, $pat_str);
    my %processed = (); 
    my $or = '';
    push @names, $name;
    push @name_desc, $pat_desc;
    foreach my $id (keys %biosqeuences){
      my $matched = $self->match_proteome_component(pattern=>\@pats,
                                                    source_type => $type,
                                                    biosequence_name => $biosqeuences{$id}{accession},
                                                    biosequence_desc => $biosqeuences{$id}{desc});
      if ($matched == 1 ){
        $entry_cnt{'redundant'}{$name}++;
        if (not defined $processed{$line}{$biosqeuences{$id}{seq}}){
					$entry_cnt{'non-redundant'}{$name}++;
					if (defined $biosqeuence_id_obs{$id} ){
						$obs{$name}++;
					}
        }
        $processed{$line}{$biosqeuences{$id}{seq}} =1;
      }
    }
  }

  #$log->error($sql);
  my @headings;
  my %heading_defs = ( Database => 'Name of database, which collectively form the reference database for this build',
                    N_Prots => 'Total number of entries in subject database',
                    N_Obs_Prots => 'Number of proteins within the subject database to which at least one observed peptide maps',
                    N_unObs_Prots => 'Number of proteins within the subject database to which none of the peptide maps',
                    Pct_Obs => 'The percentage of the subject proteome covered by one or more observed peptides' );

  for my $col_name ( qw( Database N_Prots N_Obs_Prots Pct_Obs N_unObs_Prots Description ) ) {
    push @headings, $col_name, $heading_defs{$col_name};
  }
  my $headings = $self->make_sort_headings( headings => \@headings, default => 'Database' );
  my %result;
  @{$result{table}} = ( $headings );

  my $idx = 0;
  for my $name ( @names ) {
    my $db = $name;
    my $obs  = $obs{$name} || 0;
    my $n_entry_redundant = $entry_cnt{'redundant'}{$name} || 0;
    my $n_entry_non_redundant = $entry_cnt{'non-redundant'}{$name} || 0;
    my $pct =0;
    if ( $obs && $n_entry_non_redundant ) {
        $pct = sprintf( "%0.1f", 100*($obs/$n_entry_non_redundant) );
    }
    push  @{$result{table}}, [ $db, 
                               $n_entry_redundant, 
                               $n_entry_non_redundant, 
                               $obs, 
                               $pct, 
                               $n_entry_non_redundant-$obs,
                               $name_desc[$idx]];
    $idx++;
  }
  return '' if (  @{$result{table}} == 1);
  my $table = '<table>';

  $table .= $self->encodeSectionHeader(
      text => 'Proteome Coverage (exhaustive)',
      no_toggle => 1,
      LMTABS => 1
  );

  $table .= $self->encodeSectionTable( rows => $result{table}, 
				       header => 1, 
				       table_id => 'proteome_cover',
				       align => [ qw(left left left left left left left left left left left ) ],
				       has_key => 1,
               nowrap => [qw(5 6 7 8 9 10)], 
				       rows_to_show => 25,
				       sortable => 1 );
  $table .= '</table>';

  if ($data_only){
    shift @{$result{table}};
    my $i=0;
    foreach my $row(@{$result{table}}){
      my ($org_id, $name, $type, $pat_str)  = split(/\t/, $patterns[$i]);
      $row->[0] .="|$type,$pat_str";
      $i++;
    } 
    unshift @{$result{table}} ,[qw(Database N_Entries N_Seqs N_Obs_Prots Pct_Obs N_unObs_Prots Description)]; 
    return \%result;
  }else{
    return $table;
  }
}
sub get_ptm_coverage {

  my $self = shift;
  my $build_id = shift; 
  my $patterns = shift;
  my $data_only = shift;
  my @patterns = @$patterns;

  return [] if (! $build_id || ! @patterns);
  my $sbeams = $self->getSBEAMS();
  my $sql = '';
  my $obs_sql = '';
  my @names = ();
  my %biosqeuences = ();
  my %entry_cnt = ();
  my %obs=();
  my $sql = qq~
		  SELECT B.BIOSEQUENCE_ID,B.BIOSEQUENCE_NAME, B.BIOSEQUENCE_Desc, B.biosequence_seq
			FROM $TBAT_ATLAS_BUILD AB
			JOIN $TBAT_BIOSEQUENCE B ON B.BIOSEQUENCE_SET_ID = AB.BIOSEQUENCE_SET_ID
			WHERE atlas_build_id = $build_id
      AND B.BIOSEQUENCE_NAME not like 'DECOY%'
  ~;
  my @rows = $sbeams->selectSeveralColumns($sql);
  my %biosqeuence_id_ptm =();
  foreach my $row(@rows){
    my ($id, $name, $desc,$seq) = @$row;
    $biosqeuences{$id}{accession} = $name;
    $biosqeuences{$id}{desc} = $desc;
    $biosqeuences{$id}{seq} = $seq;
  }
  $sql = qq~
    SELECT distinct ptm_type
    FROM $TBAT_PTM_SUMMARY PS
    WHERE atlas_build_id = $build_id
  ~;
  my @ptm_types = $sbeams->selectOneColumn($sql);

  my %result =();
  my $table = '';
  my %col_names=(); 
  my $counter = 0;
  foreach my $ptm_type (@ptm_types){
    $counter++;
    my %ptm_summary =();
		my $ptm_sql = qq~
				SELECT biosequence_id,
							 offset,
							 nobs,
							 residue,
							 nP95,
							 nP99,
							 nP100,      
							 nochoice, 
							 ptm_type
				FROM $TBAT_PTM_SUMMARY PS
				WHERE atlas_build_id = $build_id
				and ptm_type = '$ptm_type'
				~;
		 @rows = $sbeams->selectSeveralColumns($ptm_sql);
		 return if (! @rows);
		 my %biosqeuence_id_ptm_obs = ();
		 my %ptm_residues = ();
		foreach my $row(@rows){
			my ($id, $offset,$obs, $residue,$np95,$np99, $np100,$nochoice,$ptm_type) = @$row;
			$ptm_residues{$residue} =1;
			$biosqeuence_id_ptm_obs{$id}{$residue}{total}++; 
			if ($obs > 0){
				$biosqeuence_id_ptm_obs{$id}{$residue}{obs}++;
			}
			$biosqeuence_id_ptm_obs{$id}{$residue}{obsL1}++ if ($np95);
			$biosqeuence_id_ptm_obs{$id}{$residue}{obsL2}++ if ($np99);
			$biosqeuence_id_ptm_obs{$id}{$residue}{obsL3}++ if ($np100);
			$biosqeuence_id_ptm_obs{$id}{$residue}{obsL4}++ if ($nochoice); 
		}
		my %processed;
		foreach my $line (@patterns){
			my ($org_id, $name, $type, $pat_str, $pat_desc)  = split(/\t/, $line);
			my @pats = split(/;/, $pat_str);
			my $or = '';
			push @names, $name if ($counter == 1);
			foreach my $id (keys %biosqeuences){
				my $matched = $self->match_proteome_component(pattern=>\@pats,
																											source_type => $type,
																											biosequence_name => $biosqeuences{$id}{accession},
																											biosequence_desc => $biosqeuences{$id}{desc});
				if ($matched == 1 ){
					next if(defined $processed{$line}{$biosqeuences{$id}{seq}});
					$processed{$line}{$biosqeuences{$id}{seq}} =1;

					foreach my $site (keys %ptm_residues){
            next if ($site eq 'n');
						my @m = $biosqeuences{$id}{seq} =~ /$site/g;
						$biosqeuence_id_ptm{$id}{$site} = scalar @m;
					}

					foreach my $residue (sort {$a cmp $b} keys %ptm_residues){
						$ptm_summary{$name}{"total_$residue"} += $biosqeuence_id_ptm{$id}{$residue} || '';
					}
					if (defined $biosqeuence_id_ptm_obs{$id}){
						 foreach my $residue (keys %{$biosqeuence_id_ptm_obs{$id}}){
							 $ptm_summary{$name}{"prot_$residue"} += $biosqeuence_id_ptm_obs{$id}{$residue}{total};
							 $ptm_summary{$name}{"covered_$residue"} += $biosqeuence_id_ptm_obs{$id}{$residue}{obs};
							 $ptm_summary{$name}{"obsL1_$residue"} += $biosqeuence_id_ptm_obs{$id}{$residue}{obsL1};
							 $ptm_summary{$name}{"obsL2_$residue"} += $biosqeuence_id_ptm_obs{$id}{$residue}{obsL2};
							 $ptm_summary{$name}{"obsL3_$residue"} += $biosqeuence_id_ptm_obs{$id}{$residue}{obsL3};
							 $ptm_summary{$name}{"obsL4_$residue"} += $biosqeuence_id_ptm_obs{$id}{$residue}{obsL4};

						 }
					} 
				}
			}
		}

		#$log->error($sql);
		push @{$col_names{$ptm_type}}, 'Database';
		foreach my $residue (sort {$a cmp $b} keys %ptm_residues){
			push @{$col_names{$ptm_type}}, ("total_".$residue."_sites", "obs_protein_$residue"."_sites", "covered_$residue"."_sites"
												,"nP.80-.95","nP.95-.99","nP.99-1","no-choice");
		}
		my @headings = @{$col_names{$ptm_type}};
		my $headings = $self->make_sort_headings( headings => \@headings);
	  push @{$result{$ptm_type}}, $headings;

		for my $name ( @names ) {
			my $db = $name;
			my $obs  = $obs{$name} || 0;
			my $n_entry = $entry_cnt{$name} || 0;
			my $pct =0;
			if ( $obs && $n_entry ) {
					$pct = sprintf( "%0.1f", 100*($obs/$n_entry) );
			}
			my @summary =();
			foreach my $residue (sort {$a cmp $b} keys %ptm_residues){
				my $total = $ptm_summary{$name}{"total_$residue"} || '-'; 
				my $prot_total = $ptm_summary{$name}{"prot_$residue"}  || '-';
				my $obs = $ptm_summary{$name}{"covered_$residue"} || '-';
				if ($prot_total ne '-' ){
					push @summary, $total;
					push @summary,$prot_total;
					push @summary, sprintf("%d \(%.2f%\)",  $obs, $obs*100/$prot_total);
					push @summary,  sprintf("%d\t%d\t%d\t%d", $ptm_summary{$name}{"obsL1_$residue"},
																			$ptm_summary{$name}{"obsL2_$residue"},
																			$ptm_summary{$name}{"obsL3_$residue"},
																			$ptm_summary{$name}{"obsL4_$residue"});  
				}else{
					push @summary, ($total, '-', '-', '-','-','-','-');
				}
			}
			push @{$result{$ptm_type}}, [ $db, @summary];
		}
		#return '' if ( @result == 1);
		$table .= '<table>';

		$table .= $self->encodeSectionHeader(
				text => 'PTM Coverage',
				no_toggle => 1,
				LMTABS => 1
		);
		$table .= $self->encodeSectionTable( rows => $result{$ptm_type}, 
								 header => 1, 
								 table_id => 'ptm_coverage',
								 align => [ qw(left left left left left left left left left left left ) ],
								 has_key => 1,
								 nowrap => [qw(5 6 7 8 9 10)], 
								 rows_to_show => 25,
								 sortable => 1 );
		$table .= '</table>';
  }

  if ($data_only){
    foreach my $ptm_type (sort {$a cmp $b} keys %result){
			shift @{$result{$ptm_type}};
			my $i =0;
			foreach my $row(@{$result{$ptm_type}}){
				my ($org_id, $name, $type, $pat_str)  = split(/\t/, $patterns[$i]); 
				$row->[0] .="|$type,$pat_str"; 
				$i++;
			}
			unshift @{$result{$ptm_type}} ,[@{$col_names{$ptm_type}}];
    }
    return \%result;
  }else{
    return $table;
  }
}

sub get_proteome_coverage {
  my $self = shift;
  my $build_id = shift || return [];
  my $data_only = shift;
  my $sbeams = $self->getSBEAMS();

  my $sql = qq~
  select * 
  FROM(
		(
		SELECT COUNT(*) AS CNT, B.DBXREF_ID AS ID , DBXREF_NAME AS NAME , B.BIOSEQUENCE_SET_ID AS SETID
		FROM $TBAT_ATLAS_BUILD AB
		JOIN $TBAT_BIOSEQUENCE B ON B.biosequence_set_id = AB.biosequence_set_id
		JOIN $TBBL_DBXREF DX ON DX.dbxref_id = B.dbxref_id
		WHERE atlas_build_id = $build_id
		AND B.dbxref_id IS NOT NULL
		GROUP BY B.dbxref_id, dbxref_name, B.biosequence_set_id
		)
		UNION 
		(SELECT COUNT(DISTINCT B2.BIOSEQUENCE_ID) AS CNT, 
					 '' AS ID, 
					 LEFT (B2.BIOSEQUENCE_NAME, PATINDEX('%[_]%', B2.BIOSEQUENCE_NAME)) AS NAME, 
					 B2.BIOSEQUENCE_SET_ID AS SETID
		FROM $TBAT_ATLAS_BUILD AB2
		JOIN $TBAT_BIOSEQUENCE B2 ON B2.BIOSEQUENCE_SET_ID = AB2.BIOSEQUENCE_SET_ID
		WHERE atlas_build_id = $build_id 
		AND B2.DBXREF_ID IS NULL
		AND B2.BIOSEQUENCE_NAME LIKE '%\\_%' ESCAPE '\\'
		AND B2.BIOSEQUENCE_NAME NOT LIKE 'DECOY%'
		AND B2.BIOSEQUENCE_NAME NOT LIKE 'ENST%'
     AND B2.BIOSEQUENCE_NAME NOT LIKE 'CONTRIB%'
		GROUP BY B2.BIOSEQUENCE_SET_ID , LEFT (B2.BIOSEQUENCE_NAME, PATINDEX( '%[_]%', B2.BIOSEQUENCE_NAME ))
		)
    UNION
    (SELECT COUNT(DISTINCT B2.BIOSEQUENCE_ID) AS CNT,
           '' AS ID,
           LEFT (B2.BIOSEQUENCE_NAME, charindex( '_' , B2.BIOSEQUENCE_NAME,charindex ('_',B2.BIOSEQUENCE_NAME, 0) +1) ),
           B2.BIOSEQUENCE_SET_ID AS SETID
    FROM $TBAT_ATLAS_BUILD AB2
    JOIN $TBAT_BIOSEQUENCE B2 ON B2.BIOSEQUENCE_SET_ID = AB2.BIOSEQUENCE_SET_ID
    WHERE atlas_build_id = $build_id
    AND B2.DBXREF_ID IS NULL
     AND B2.BIOSEQUENCE_NAME LIKE 'CONTRIB%'
    AND B2.BIOSEQUENCE_NAME NOT LIKE 'CONTRIB_ENST%'
    GROUP BY B2.BIOSEQUENCE_SET_ID ,LEFT (B2.BIOSEQUENCE_NAME, charindex( '_' , B2.BIOSEQUENCE_NAME,charindex ('_',B2.BIOSEQUENCE_NAME, 0) +1) ) 
    )
    UNION
    (SELECT COUNT(DISTINCT B2.BIOSEQUENCE_ID) AS CNT,
           '' AS ID,
           'CONTRIB_ENST_' AS NAME, 
           B2.BIOSEQUENCE_SET_ID AS SETID
    FROM $TBAT_ATLAS_BUILD AB2
    JOIN $TBAT_BIOSEQUENCE B2 ON B2.BIOSEQUENCE_SET_ID = AB2.BIOSEQUENCE_SET_ID
    WHERE atlas_build_id = $build_id
    AND B2.DBXREF_ID IS NULL
    AND B2.BIOSEQUENCE_NAME LIKE 'CONTRIB_ENST%'
    GROUP BY B2.BIOSEQUENCE_SET_ID
    )
  )CAT
  ORDER BY NAME
  ~;
  my $sth = $sbeams->get_statement_handle( $sql );
  my @names;
  while ( my @row = $sth->fetchrow_array() ) {
    push @names, \@row;
  }

  my $obs_sql = qq~
  (SELECT COUNT(DISTINCT B.BIOSEQUENCE_ID),convert(varchar(10), B.DBXREF_ID)
  FROM $TBAT_BIOSEQUENCE B
  JOIN $TBAT_PEPTIDE_MAPPING PM 
    ON PM.MATCHED_BIOSEQUENCE_ID = B.BIOSEQUENCE_ID
  JOIN $TBAT_PEPTIDE_INSTANCE PI 
    ON PM.PEPTIDE_INSTANCE_ID = PI.PEPTIDE_INSTANCE_ID
  WHERE atlas_build_id = $build_id
  GROUP BY B.DBXREF_ID
  )
  UNION
  (SELECT COUNT(DISTINCT B2.BIOSEQUENCE_ID) AS CNT, 
         LEFT (B2.BIOSEQUENCE_NAME, PATINDEX('%[_]%', B2.BIOSEQUENCE_NAME)) AS CAT
    FROM $TBAT_BIOSEQUENCE B2
    JOIN $TBAT_PEPTIDE_MAPPING PM2
      ON PM2.MATCHED_BIOSEQUENCE_ID = B2.BIOSEQUENCE_ID
    JOIN $TBAT_PEPTIDE_INSTANCE PI2
      ON PM2.PEPTIDE_INSTANCE_ID = PI2.PEPTIDE_INSTANCE_ID
    WHERE atlas_build_id = $build_id
  and B2.dbxref_id is null
  AND B2.BIOSEQUENCE_NAME LIKE '%\\_%' ESCAPE '\\'
  AND B2.BIOSEQUENCE_NAME NOT LIKE 'DECOY%'
  AND B2.BIOSEQUENCE_NAME NOT LIKE 'CONTRIB%'
  GROUP BY LEFT (B2.BIOSEQUENCE_NAME, PATINDEX( '%[_]%', B2.BIOSEQUENCE_NAME ))
  )
  UNION
  (SELECT COUNT(DISTINCT B2.BIOSEQUENCE_ID) AS CNT,
         LEFT (B2.BIOSEQUENCE_NAME, charindex( '_' , B2.BIOSEQUENCE_NAME,charindex ('_',B2.BIOSEQUENCE_NAME, 0) +1) ) AS CAT
    FROM $TBAT_BIOSEQUENCE B2
    JOIN $TBAT_PEPTIDE_MAPPING PM2
      ON PM2.MATCHED_BIOSEQUENCE_ID = B2.BIOSEQUENCE_ID
    JOIN $TBAT_PEPTIDE_INSTANCE PI2
      ON PM2.PEPTIDE_INSTANCE_ID = PI2.PEPTIDE_INSTANCE_ID
    WHERE atlas_build_id = $build_id
  and B2.dbxref_id is null
  AND B2.BIOSEQUENCE_NAME LIKE 'CONTRIB%'
  AND B2.BIOSEQUENCE_NAME NOT LIKE 'CONTRIB_ENST%'
  AND B2.BIOSEQUENCE_NAME NOT LIKE 'CONTRIB_H%'
  GROUP BY LEFT (B2.BIOSEQUENCE_NAME, charindex( '_' , B2.BIOSEQUENCE_NAME,charindex ('_',B2.BIOSEQUENCE_NAME, 0) +1) ) 
  )
  UNION
  (SELECT COUNT(DISTINCT B2.BIOSEQUENCE_ID) AS CNT,
         'CONTRIB_ENST_' AS CAT 
    FROM $TBAT_BIOSEQUENCE B2
    JOIN $TBAT_PEPTIDE_MAPPING PM2
      ON PM2.MATCHED_BIOSEQUENCE_ID = B2.BIOSEQUENCE_ID
    JOIN $TBAT_PEPTIDE_INSTANCE PI2
      ON PM2.PEPTIDE_INSTANCE_ID = PI2.PEPTIDE_INSTANCE_ID
    WHERE atlas_build_id = $build_id
  and B2.dbxref_id is null
  AND B2.BIOSEQUENCE_NAME LIKE 'CONTRIB_ENST%'
  )

  ~;

  my $obssth = $sbeams->get_statement_handle( $obs_sql );
  my %obs;
  while ( my @row = $obssth->fetchrow_array() ) {
    $obs{$row[1]} = $row[0];
  }

  my @headings;
  my %head_defs = ( Database => 'Name of databse, which collectively form the reference database for this build',
                    N_Seqs => 'Total number of entries in subject database',
                    N_Obs_Prots => 'Number of proteins within the subject database to which at least one observed peptide maps',
                    Pct_Obs => 'The percentage of the subject proteome covered by one or more observed peptides' );

  for my $head ( qw( Database N_Seqs N_Obs_Prots Pct_Obs ) ) {
    push @headings, $head, $head_defs{$head};
  }
  my $headings = $self->make_sort_headings( headings => \@headings, default => 'Database' );
  my @return = ( $headings );

  for my $row ( @names ) {
    my $db = $row->[2]; 
    $db =~ s/\_$//;
    if ( $db eq 'Swiss-Prot' ) {
      $db .= ' (may include Varsplic)';
    } elsif ( $db eq 'UniProt' ) {
      $db .= ' (excludes SwissProt if shown)';
    }
    my $pct = 0;
    my $obs = '';
    if($row->[1] =~ /^\d+$/ && $row->[1] !~ /^0$/){ 
      $obs = $obs{$row->[1]};
    }else{
      $obs = $obs{$row->[2]};
    }

    if ( $obs && $row->[0] ) {
        $pct = sprintf( "%0.1f", 100*($obs/$row->[0]) );
    }
    push @return, [ $db, $row->[0], $obs, $pct ];
  }
  return '' if ( @return == 1);
  my $table = '<table>';

  $table .= $self->encodeSectionHeader(
      text => 'Proteome Coverage (exhaustive)',
      no_toggle => 1,
      LMTABS => 1
  );

  $table .= $self->encodeSectionTable( rows => \@return, 
				       header => 1, 
				       table_id => 'proteome_cover',
				       align => [ qw(left right right right ) ], 
				       rows_to_show => 25,
				       has_key => 1,
				       sortable => 1 );
  $table .= '</table>';

  if ($data_only){
    return \@return;
  }else{
    return $table;
  }
}

sub get_what_is_new {
  my $self = shift;
  my $build_id = shift || return [];
  my $data_only = shift;

  my $sbeams = $self->getSBEAMS();
  # check if it is default build 
  my $sql = qq~
    SELECT  DEFAULT_ATLAS_BUILD_ID
    FROM $TBAT_DEFAULT_ATLAS_BUILD 
    WHERE ATLAS_BUILD_ID = $build_id
  ~;
  my @row = $sbeams->selectOneColumn($sql); 
  return if(! @row);
    
  # compare biosequence set version
  $sql = qq~
    SELECT AB.ATLAS_BUILD_NAME, BS.SET_DESCRIPTION    
    FROM $TBAT_ATLAS_BUILD AB
    JOIN $TBAT_BIOSEQUENCE_SET BS ON (BS.BIOSEQUENCE_SET_ID = AB.BIOSEQUENCE_SET_ID ) 
    WHERE AB.ATLAS_BUILD_ID = $build_id     
  ~;
  @row = $sbeams->selectSeveralColumns($sql);
  my ($default_build_name, $default_bsset_desc) = @{$row[0]};
  my $build_name_pat = $default_build_name;
  $build_name_pat =~ s/\s+\d.*//;
  $sql = qq~
    SELECT AB.atlas_build_id, AB.ATLAS_BUILD_NAME, BS.SET_DESCRIPTION
    FROM $TBAT_ATLAS_BUILD AB
    JOIN $TBAT_BIOSEQUENCE_SET BS ON (BS.BIOSEQUENCE_SET_ID = AB.BIOSEQUENCE_SET_ID )
    WHERE AB.ATLAS_BUILD_NAME like '$build_name_pat [0-9]%' and AB.project_id = 475 
    AND AB.ATLAS_BUILD_ID != $build_id
  ~;
  @row = $sbeams->selectSeveralColumns($sql);
  
  return if (! @row);
  my ($previous_build_id, $previous_build_name, $previous_bsset_desc) = @{$row[0]};

  my @headings;
  push @headings, '', '';
  push @headings, $default_build_name, '';
  push @headings ,  $previous_build_name, '';

  my $headings = $self->make_sort_headings( headings => \@headings, default => $default_build_name );
  my @return = ( $headings );
  push @return , ['Reference_Database', $default_bsset_desc, $previous_bsset_desc];

  ## get sample number
  $sql = qq~ 
      SELECT ATLAS_BUILD_ID , COUNT(SAMPLE_ID)  
      FROM $TBAT_ATLAS_BUILD_SAMPLE 
      WHERE atlas_build_id in ($build_id, $previous_build_id)
      GROUP BY ATLAS_BUILD_ID
  ~;
  my %sample_cnt = $sbeams->selectTwoColumnHash($sql);
  push @return , ['# Experiments'];
  push @return , ['Distinct_Peptides'];
  push @return , ['Canonical_Proteins'];

  foreach my $bid (($build_id, $previous_build_id)){
    $sbeams->update_PA_table_variables($bid); 
		$sql = qq~
		SELECT  atlas_build_id,  COUNT(distinct PI.peptide_instance_id) cnt
		FROM $TBAT_PEPTIDE_INSTANCE PI
		JOIN $TBAT_PEPTIDE_MAPPING PM ON (PI.PEPTIDE_INSTANCE_ID = PM.PEPTIDE_INSTANCE_ID)
		JOIN $TBAT_BIOSEQUENCE B ON (PM.MATCHED_BIOSEQUENCE_ID = B.BIOSEQUENCE_ID)
		WHERE atlas_build_id in ($bid)
		AND B.BIOSEQUENCE_NAME NOT LIKE 'DECOY%'
		AND B.BIOSEQUENCE_NAME NOT LIKE 'CONTAM%'
		AND B.biosequence_name NOT LIKE '%UNMAPPED%'
		AND B.biosequence_desc NOT LIKE '%common contaminant%'
		GROUP BY atlas_build_id
		~;
		my %pep_count = $sbeams->selectTwoColumnHash($sql);

		$sql = qq~ 
		SELECT PID.atlas_build_id, COUNT(BS.biosequence_name) cnt
		FROM $TBAT_PROTEIN_IDENTIFICATION PID
		JOIN $TBAT_PROTEIN_PRESENCE_LEVEL PPL
		ON PPL.protein_presence_level_id = PID.presence_level_id
		JOIN $TBAT_BIOSEQUENCE BS
		ON BS.biosequence_id = PID.biosequence_id
		WHERE PID.atlas_build_id in ($bid)
		AND PPL.level_name = 'canonical'
		AND BS.biosequence_name NOT LIKE 'DECOY%'
		AND BS.biosequence_name NOT LIKE '%UNMAPPED%'
		AND BS.biosequence_desc NOT LIKE '%common contaminant%'
		AND BS.BIOSEQUENCE_NAME NOT LIKE 'CONTAM%'
		GROUP BY PID.atlas_build_id
		~;
		my %prot_count = $sbeams->selectTwoColumnHash($sql);

		push @{$return[2]}, $sample_cnt{$bid};
		push @{$return[3]}, $pep_count{$bid};
		push @{$return[4]}, $prot_count{$bid};
  }
  my $table = '<table>';
  $table .= $self->encodeSectionHeader(
      LMTABS => 1,
      no_toggle => 1,
      text => "What&#39s New",
  );

  $table .= $self->encodeSectionTable( rows => \@return,
				       header => 1,
				       table_id => 'what_is_new',
				       align => [ qw(left right right right ) ],
				       bg_color => '#f3f1e4', #EAEAEA',
				       has_key => 1,
				       rows_to_show => 25,
				       sortable => 1 );
  $table .= '</table>';
  ## new sample table:
  print "current_build=$build_id\n";
  print "previous_build_id=$previous_build_id\n";
  $sql = qq~
    SELECT SAMPLE_ID
    FROM $TBAT_SAMPLE 
    WHERE REPOSITORY_IDENTIFIERS IN ( 
      SELECT SA.REPOSITORY_IDENTIFIERS 
      FROM $TBAT_ATLAS_BUILD_SAMPLE A 
      JOIN $TBAT_SAMPLE SA ON (A.sample_id = SA.sample_id)
      WHERE A.ATLAS_BUILD_ID IN ($build_id)
    )
    AND SAMPLE_ID in (
     SELECT A.sample_id 
      FROM $TBAT_ATLAS_BUILD_SAMPLE A
      WHERE A.ATLAS_BUILD_ID IN ($build_id)
    )
    AND REPOSITORY_IDENTIFIERS NOT IN ( 
      SELECT SB.REPOSITORY_IDENTIFIERS 
      FROM $TBAT_ATLAS_BUILD_SAMPLE B
      JOIN $TBAT_SAMPLE SB ON (B.sample_id = SB.sample_id) 
      WHERE B.ATLAS_BUILD_ID IN ($previous_build_id)
    )
    AND SAMPLE_ID not IN(
     SELECT B.sample_id
      FROM $TBAT_ATLAS_BUILD_SAMPLE B
      WHERE B.ATLAS_BUILD_ID IN ($previous_build_id)
    )
    order by REPOSITORY_IDENTIFIERS 
  ~;

  my @sample_ids = $sbeams->selectOneColumn($sql);
  if(@sample_ids){
    my $sampleDisplay = $self->getProteinSampleDisplay( sample_ids => \@sample_ids,
																											no_header => 1,
																											rows_to_show => 25,
																											sample_category => 1,
																											#bg_color  =>  '#f3f1e4', #EAEAEA'
																											#sortable => 1,
																											max_rows => 500);
    $table .=$sbeams->make_toggle_section( neutraltext => 'New Experiments',
					   sticky => 1,
					   name => 'getnew_samplelist_div',
					   barlink => 1,
					   visible => 1,
					   content => "<table>$sampleDisplay</table>" );
  }

  $sbeams->update_PA_table_variables($build_id);
  if ($data_only){
    return (\@return,  \@sample_ids);
  }else{ 
    return $table;
  }
}

sub get_scroll_table {
  my $self = shift || die ("self not passed");
  my %args = ( width => 900,
               bold => 1,
               key_width => 20,
               @_ );

  for my $arg ( qw( sql headings ) ) {
    die "Missing required argument $arg" if !defined $args{$arg};
  }

  $sbeams = $self->getSBEAMS();
  my $sth = $sbeams->get_statement_handle->( $args{sql} );
  while ( my @row = $sth->fetchrow_array() ) {
  }

}

sub print_html_table_to_tsv{
  my $self = shift;
  my %args = @_;
  my $data_ref = $args{data_ref};
  my $column_name_ref = $args{column_name_ref};
  my $filename = $args{filename};
  my %params =();
  $params{'Content-Disposition'}="attachment;filename=$args{filename}";
  print $q->header(-type=>'tsv',%params);

  foreach (@$column_name_ref){
    print lc($_) ."\t";
  }
  print "\n";
  foreach my $row (@$data_ref){
    foreach my $val (@$row){
      $val =~ s/(\<[^\<]+\>)/ /g;
      print "$val\t";
    }
    print "\n";
  }

}
sub display_chromosome_coverage_plotly{
  my $self = shift;
  my %args = @_;
  my $data_ref = $args{data_ref};
  ## column names
  #1. organism
  #2. chromosome
  #3. total count
  #4. observed count
  my %data = ();
  my $max = 0;
  my $charts = '<div style="width: 60vw">';
  my %chromosome = ();
  foreach my $row (@$data_ref){
    my ($org, $chr,$total, $observed) = @$row;
    $chromosome{$chr} =1;
    $data{$org}{total}{$chr} = $total;
    $data{$org}{obs}{$chr} = $observed;
    #push @{$data{$org}{total}},[$chr, $total];
    #push @{$data{$org}{obs}} , [$chr,$observed];
    $max = $total if ($total > $max);
  }
  $max = sprintf("%.0f", log($max)/log(10));
  my $n_orgs = scalar keys %data;
  foreach my $org (keys %data){
    my $plot_title = '';
    my $chart;
    my (@bdata_a, @bdata_b);
    foreach my $chr (sort {$a cmp $b} keys %chromosome){
      if ($data{$org}{total}{$chr}){
        push @bdata_a, [$chr, $data{$org}{total}{$chr}];
        push @bdata_b, [$chr, $data{$org}{obs}{$chr}];
      }else{
        push @bdata_a, [$chr, ''];
        push @bdata_b, [$chr, ''];
      }
    }

    if ($n_orgs > 1){
      $chart = $self->plotly_barchart (data => [\@bdata_a, \@bdata_b],
                                       names => ['total', 'observed'],
																			 divName => "chr_cov_div_$org",
																			 ytitle => 'Count',
                                       dtick => 1,
                                       barlabel => 1,
                                       xtickangle => 'tickangle:25',
                                       yticktype => "type:'log',range:[0,$max],dtick:1",
                                       layoutmargin => 'b:50',
                                       title => $org);
     }else{
       $chart = $self->plotly_barchart (data => [\@bdata_a, \@bdata_b],
                                       names => ['total', 'observed'],
                                       divName => "chr_cov_div_$org",
                                       ytitle => 'Count',
                                       dtick => 1,
                                       xtickangle => 'tickangle:25',
                                       yticktype => "type:'log',range:[0,$max],dtick:1",
                                       layoutmargin => 'b:50',
                                       barlabel => 1
                                       );
     }
    $charts .= qq~
      <div id="chr_cov_div_$org"</div>
      $chart
    ~;   
  }

  my  $html .= $sbeams->make_toggle_section(
      neutraltext =>"Genome Coverage",
      sticky => 1,
      barlink => 1,
      visible => 1,
      name => "chr_plots_div",
      content => "<TABLE><TR><TD COLSPAN='5'>$charts</div></TD></TR></TABLE>",
  );
 
  return "$html";

}


sub display_peptide_sample_category_plotly{
  my $self = shift;
  my %args = @_;
  my $data_ref = $args{data_ref};
  my $sample_array_ref = $args{sample_array_ref};
  my $build_id = $args{build_id};
  my $column_name_ref = $args{column_name_ref};

  my %cols=();
  my $idx = 0;
  foreach (@$column_name_ref){
    $cols{$_} = $idx;
    $idx++;
  }
  my (@sample_category, @peptide_count,@links, @pect_spec_per_sample );
  my $total_observed_spectra;
  my @n_good_spectra; 

 foreach my $data (@$data_ref){
    my ($name,$id,$cnt) = @$data;
    my ($good_spectra) =0;
    foreach my $row (@$sample_array_ref){
      my $n = scalar @$row; 
      my $sample_cat_id  = $row->[$cols{"sample_category_id"}];
      $row->[$cols{"Spectra ID'd"}] =~ s/,//g;
      if ($sample_cat_id eq $id){
        $good_spectra +=  $row->[$cols{"Spectra ID'd"}];
        $total_observed_spectra += $row->[$cols{"Spectra ID'd"}]; 
      }
    }  
    push @n_good_spectra, $good_spectra;
    push @peptide_count, $cnt; 
    push @sample_category,'<a href="'."https://db.systemsbiology.net/sbeams/cgi/PeptideAtlas/GetPeptides?atlas_build_id=$build_id&sample_category_id=$id&QUERY_NAME=AT_GetPeptides&apply_action=QUERY" .'">'. $name .'  </a>';
  }
 
  foreach my $val (@n_good_spectra){
    push @pect_spec_per_sample,  sprintf( "%0.6f", ($val/$total_observed_spectra)*100 );
  } 

  my $sample_category_str = join("','", @sample_category);
  my $n_good_spectra_str = join(",", @n_good_spectra);  
  my $pect_spec_per_sample_str = join(",", @pect_spec_per_sample);
  my $nc = $#sample_category + 1;
  my $height = $nc*40 > 1100? 1100:$nc*40; 
  $height = $height < 300?300:$height;
  my $plot_js = qq~
			var pepcnt = {
				y: ['$sample_category_str'],
				x: [$pect_spec_per_sample_str],
        //customdata:['link_str'],
				name:'Percentage of PSMs per sample category',
			  type: 'bar',
        marker: {color: '0x5588bb'},	
        orientation: 'h',
        opacity: 0.8
			};
      var obs = {
        y: ['$sample_category_str'],
        x: [$n_good_spectra_str],
        name:'Total PSMs per sample category',
        type: 'bar',
        marker: {color: '0x5588bb'},
        orientation: 'h',
        opacity: 0.8,
      };

			var layout = {
        legend:{x: 0.029,y: 1.1},
			  height: $height,
        font: {size: 12},
				title:'Percentage of PSMs per sample category',
        margin:{l: 350},
        hoverlabel:{bgcolor:'white'},
			};
			var data1 = [pepcnt];
      var config = {
        toImageButtonOptions: {
          format: 'png', // one of png, svg, jpeg, webp
          filename: 'image',
          height: 800,
          width: 1000,
          scale: 6,
        }
      };

			Plotly.newPlot('plot_div3', data1,layout,config);
      layout.title = "Total PSMs per sample category";
      var data2 = [obs];
      Plotly.newPlot('plot_div4', data2,layout,config);

  ~;
  my $chart = qq~
    <script type="text/javascript" src="https://cdn.plot.ly/plotly-latest.min.js"></script>
	  <script type='text/javascript'>
			function toggle_plot() {
        var myplot3 =  document.getElementById('plot_div3');
        var myplot4 =  document.getElementById('plot_div4');
        if (myplot3.style.display == '' ){
           myplot3.style.display = 'none';
           myplot4.style.display = '';
          document.getElementById('toggle_button').innerHTML = 'Show Percentage of PSMs per sample category plot';
        } else {
           myplot3.style.display = '';
           myplot4.style.display = 'none';
          document.getElementById('toggle_button').innerHTML = 'Show Total PSMs per sample category plot';
        }
			}
		</script>
       <A NAME='<B>Sample Category</B>'></A><DIV CLASS="hoverabletitle"i onclick="toggle_content('sample_cat_div')">
       <IMG ID='sample_cat_div_gif'  SRC='/devZS/sbeams/images/small_gray_minus.gif'>
       <B>Spectra Identification by Sample Category</B></DIV>
       <DIV ID='sample_cat_div' class="visible">
       <TABLE><TR><TD COLSPAN='5'>
       <TR> 
        <TD><button type='button' id='toggle_button' onclick=toggle_plot()>Show Percentage of PSMs per sample category plot</button>
           <div style="width: 80vw">
            &nbsp;<div id="plot_div3" style="display:none; width: 80vw"> </div>
            &nbsp;<div id="plot_div4"></div>
            <br><br>
           <div>
       </TD>
       </TR>
       </TD></TR></TABLE></DIV> 
		<script type="text/javascript" charset="utf-8">
      $plot_js
		</script>
  ~;
  return $chart;
  
}
sub displayProt_PTM_plotly{
   my $self = shift;
  my %args = @_;
  my $protein = $args{protein};
  my $data = $args{data};
  my $sequence = $args{seq},
  my $atlas_build_id = $args{atlas_build_id};
  my @column_names = @{$args{column_names}};
  my $ptm_obs = $args{ptm_obs}; 
  my @ptm_types = keys %$data; 
  my $div_counter=1;
  my $plot_js =qq~
      function getCol(matrix, col){
         var column = [];
         for(var i=0; i<matrix.length; i++){
            column.push(matrix[i][col]);
         }
         return column;
      };
      var myPlots=[];
      var allData=[];
      var allObs=[];
      var allXvals=[];
      var allTicklabels=[];
   ~;

	my $seqlen = length($sequence);
  foreach my $ptm_type (@ptm_types){
		my $dataTable = "var data$div_counter=[['','','','','','','','',''],\n";
		my $totalObs = "var totalObs$div_counter=[";
		my $maxnobs  = 0;
		my $tickvals   = "var tickvals$div_counter =['0'"; 
		#$sequence =~ s/\*.*//g;
		my @aas = split(//, $sequence);
		my $sep = ',';
		foreach my $pos (0..$#aas){
			my $aa=$aas[$pos];
			my $nobs = '';
			#if ($aa =~ /[$residue]/i){
			if (defined $data->{$ptm_type}{$protein}{$pos}){
        my $flag = 0;
        foreach my $row(@{$data->{$ptm_type}{$protein}{$pos}}){ 
          $nobs ='';
          $flag = 1 if($row->[0] !~ /[nc]/);
					$dataTable .= "['$row->[0]'";
					#foreach my $c (qw (nP<.01 nP<.05 nP<.20 nP.2-.8 nP>.80 nP>.95 nP>.99 no-choice)){
          foreach my $i (6..13){
						$dataTable .=",$row->[$i]";
						if ($row->[$i] > $maxnobs){
							$maxnobs =$row->[$i];
						}
						$nobs += $row->[$i]; 
					}
					$dataTable .= "],\n";
          $tickvals .= "$sep'";
					$tickvals .= $pos+1;
          $tickvals .= "<br>$row->[0]'";
					$totalObs .= "$sep";
					$totalObs .= $nobs;
        }
        ### n-terminal AA not modified 
        if (! $flag){
            $dataTable .= "['$aa','','','','','','','',''],\n";
						$tickvals .= "$sep'";
						$tickvals .= $pos+1;
            $tickvals .= "<br>$aa'";
						$totalObs .= "$sep". '';
        } 
        next;
			}
 
			$dataTable .= "['$aa','','','','','','','',''],\n";
			$tickvals .= "$sep'"; 
			$tickvals .=$pos+1;
      $tickvals .= "<br>$aa'";
			$totalObs .= "$sep";
			$totalObs .= $nobs;
		}
		$dataTable =~ s/,$//;
		$dataTable =~ s/\n$//;
		$dataTable .= "];\n";
		$tickvals   .= "];\n"; 
		$totalObs .= "];\n";

		my $trace = "var col1 =  getCol(data$div_counter, 0);\n";
		my @colors = ('red','orange','purple','grey','#007eca','skyblue','green','black');
		my @names = ('<0.01','0.01-0.05','0.05-0.20','0.20-0.80','0.80-0.95','0.95-0.99','0.99-1.00','no-choice');
		foreach my $i(1..8){
			my $j=$i+1;
			my $k=$i-1;
			$trace .= qq~
			var col$j = getCol(data$div_counter, $i);
			var trace$i = {
				x:xvals$div_counter,
				y:col$j,
				marker: {color:'$colors[$k]'},
				type: 'bar',
				name:'$names[$k]'
			};
			~;
		}

	 $plot_js .= qq~
			$dataTable;
			$tickvals
			$totalObs;
			var myDiv$div_counter = document.getElementById("protPTM_div$div_counter");
      myPlots.push (myDiv$div_counter);
      var xvals$div_counter = [];
      for(var i=0; i<tickvals$div_counter.length; i++){
        xvals$div_counter.push(tickvals$div_counter\[i]);
      }

			$trace
			var layout$div_counter = {
				annotations: [],
				hovermode:'x unified',
				barmode:'stack',
				yaxis:{title:'nObs',rangemode:'tozero'},
        xaxis:{tickvals: xvals$div_counter,
               type: 'category',
               dtick:50,
               tick0: '0',
               tickmode:'linear',
               ticktext:tickvals$div_counter},
			};
		 
			for ( var i = 0 ; i < tickvals$div_counter.length; i++ ) {
				if (totalObs$div_counter\[i] !== null && totalObs$div_counter\[i] > 0){
					var anno = {
						x: i, 
						y: totalObs$div_counter\[i], 
						text: '<b>'+totalObs$div_counter\[i]+'</b>',
						xanchor: 'center',
						yanchor: 'top',
						yshift:30,
						showarrow: false
					};
					layout$div_counter.annotations.push(anno);
			 }
		 }
		 plotdata$div_counter= [trace1,trace2,trace3,trace4,trace5,trace6,trace7,trace8];
		 Plotly.newPlot(myDiv$div_counter, plotdata$div_counter, layout$div_counter, {scrollZoom: true});
     ~;
     $div_counter++;
   }

   my $counter = 1;
   while ($div_counter > $counter){
     $plot_js .= "allData.push(data$counter);\n";
     $plot_js .= "allObs.push(totalObs$counter);\n";
     $plot_js .= "allXvals.push(xvals$counter);\n";
     $counter++;
   }

   $plot_js .=  qq~
	 function relayout(eventdata, div) {
		 //alert(JSON.stringify(eventdata));
     if (Object.entries(eventdata).length === 0) {return;}
			var xmax=eventdata['xaxis.range[1]'];
			var xmin=eventdata['xaxis.range[0]'];
			var data;
			var totalObs;
      var xvals;
      for (var i=0; i<myPlots.length; i++){
				if (div == myPlots[i] ){
						data = allData[i];
						totalObs=allObs[i];
            xvals=allXvals[i];
				}
		  }

			var n = xmin  < 0 ? 0:parseInt(xmin);
			var m = xmax  > $seqlen ? $seqlen-1 :parseInt(xmax);
			var ymax=0;
			var len = data.length;
			//get max obs
			var l =0;
			for (var i =n; i <=m; i++){
				var o = parseInt(totalObs[i]) || 0;
				l += o;
				if(o > ymax){	
					ymax=o;
				}    
			}
			ymax = ymax + 5;
			if ( xmax- xmin <= 50){
			 var update = {
				 xaxis:{
          type: 'category',
					tickangle:0,
          tickmode:'array',
					range: [xmin, xmax], 
				 },
				 yaxis:{
					title:'nObs',
					range:[0,ymax]
				 }
			 };
				Plotly.relayout(div, update);
			}
			if ( xmax- xmin > 50){	
				 if (xmin < 0){
					 xmin=0;
				 }
				 if (xmax> $seqlen){
					 xmax=$seqlen;
				 }
				 var update = {
					 xaxis:{
					  tickmode:'linear',
					  dtick:50,
            tick0: xvals[xmin],
            type: 'category',
						range: [xmin, xmax], 
					 },
					 yaxis:{
						title:'nObs',
						range:[0,ymax]
					 }
				 };
				 Plotly.relayout(div, update);
			 }
       if (eventdata['xaxis.autorange'] && eventdata['yaxis.autorange'] ){
         var update = {
           yaxis:{title:'nObs',rangemode:'tozero'},
           xaxis:{
             tickmode:'linear',
             ticktext:xvals,
             dtick:50,
             tick0:'0'
           }
         };
         Plotly.relayout(div, update);
       }

   };
	 myPlots.forEach(div => {
			div.on("plotly_relayout", function(ed) {
				relayout(ed, div);
			});
	 });
  ~;
  my $chart = qq~
    <div class="tab">
			<style>
			/* Style the tab */
			.tab {
				overflow: hidden;
			}
			/* Style the buttons inside the tab */
			.tab button {
			  padding: 14px 16px;
        border: none;
        outline: none;	
        font-size: 14px;
        font-weight:bold;
				float:left;
				color:#555555;
				text-decoration:none;
				background:#f3f1e4;
        border-right:1.25px solid #AAAAAA;       
			}

			/* Change background color of buttons on hover */
			.tab button:hover {
				background-color:#bb0000;
            color:#ffffff; 
			}

			/* Create an active/current tablink class */
			.tab button.active {
				 background:#ffffff;
			    color:#333333; 
			}

			/* Style the tab content */
			.tabcontent {
				border-top: none;
			}
			</style>
  ~;
  $counter = 1;
  foreach my $ptm_type(@ptm_types){
    $chart .="<button type='button' class='ptmtablinks' onclick='openPTM(event,\"$ptm_type\")'>$ptm_type</button>\n";
    $counter++;
  }
  $chart .="</div>";

  $counter = 1;
  foreach my $ptm_type(@ptm_types){ 
    my @names = @column_names;
    my $PTM_table_Display = $self->getPTMTableDisplay( cols => \@names,
                    ptm_type => $ptm_type,
                    data => $data->{$ptm_type},
                    atlas_build_id => $atlas_build_id,
                    biosequence_name => $protein,);
                    #rows_to_show => 10,
                    #max_rows => 500);
    @names = @column_names[1..$#column_names-1];
    my $ptm_table_help = $self->get_table_help(column_titles_ref=>\@names);
		my $text = $ptm_type;
		if ($ptm_type=~ /[STY]/i){
			 $text = "Phospho Summary ($ptm_type)" ;
		}elsif($ptm_type =~ /[nk]\:42/i){
			 $text = "Acetylation Summary ($ptm_type)";
		}elsif ($ptm_type =~ /[nKR]\:14.015/i){
       $text = "Methylation Summary ($ptm_type)";
    }elsif ($ptm_type eq ''){
      $text = "PTM Summary";
    }
    my $html = '';
    my $checked = '';
    my $hidden='true';
    if ($counter==1){
      $checked = "checked";
      $hidden='false';
    }
    if ($counter == 1){
      $chart .= "<div  id ='$ptm_type' style='display:block' class='tabcontent'>\n";
    }else{
      $chart .= "<div  id ='$ptm_type' style='display:none' class='tabcontent'>\n";
    }
    if ($ptm_obs->{$ptm_type}){
      $chart .= qq~<div style="width: 80vw" id="protPTM_div$counter"></div>\n~;
    }else{
      $chart .= qq~<div style="width: 80vw;display:none" id="protPTM_div$counter"></div>\n~;
    }
    $chart .= qq~$ptm_table_help\n<TABLE width='100%'>$PTM_table_Display</TABLE>\n~;
     $counter++;
    $chart .= "</div>";

  }
	$chart = $sbeams->make_toggle_section( neutraltext => "PTM Summary",
	             sticky => 1,
	             name => "getprotein_ptm_summary_div",
	             barlink => 1,
	             visible => 1,
	             content => $chart);
  $chart .=  qq~
     <!-- Latest compiled and minified plotly.js JavaScript -->
     <script type="text/javascript" src="https://cdn.plot.ly/plotly-latest.min.js"></script>
     <script type="text/javascript" charset="utf-8">
			function openPTM(evt, Name) {
				var i, tabcontent, ptmtablinks;
				tabcontent = document.getElementsByClassName("tabcontent");
				for (i = 0; i < tabcontent.length; i++) {
					tabcontent[i].style.display = "none";
				}
				ptmtablinks = document.getElementsByClassName("ptmtablinks");
				for (i = 0; i < ptmtablinks.length; i++) {
					ptmtablinks[i].className = ptmtablinks[i].className.replace(" active", "");
				}
				document.getElementById(Name).style.display = "block";
				evt.currentTarget.className += " active";
			}
      $plot_js
     </script>
  ~;
  
  return $chart;

}
sub displayExperiment_contri_plotly{
  my $self = shift;
  my %args = @_;
  my $data_ref = $args{data_ref};
  my $tr = $args{tr};
  my $column_name_ref = $args{column_name_ref};
  my %cols =();
  my @sample_label = ();
  my (@cumpepx, @cumpepy,@idvpepy,@cumprotx, @cumproty,@idvproty);
  my $pre_cum_n_good_spectra;
  my $idx =0;
  foreach (@$column_name_ref){
    $cols{$_} = $idx;
    $idx++;
  } 
  $idx=0;
  
  foreach my $row(@$data_ref){
    my $n_good_spectra = $row->[$cols{"Spectra ID'd"}];
    my $n_distinct_peptides = $row->[$cols{"Distinct Peptides"}];
    my $cumulative_n_peptides = $row->[$cols{"Cumulative Peptides"}];
    my $n_canonical_proteins = $row->[$cols{"Distinct Canonical Proteins"}];
    my $cumulative_n_proteins = $row->[$cols{"Cumulative Canonical Proteins"}];
    my $sample_tag =  $row->[$cols{"Experiment Tag"}];
    $sample_tag =~ s/[\n\r]//g;
    $sample_tag =~ s/.*sample_id=\d+["']\s?>//;
    $sample_tag =~ s/<.*//;
    push @sample_label, ($sample_tag,'','');
    $n_good_spectra =~ s/,//g;
    $n_distinct_peptides =~ s/,//g;
    $cumulative_n_peptides  =~ s/,//g;
    $n_canonical_proteins  =~ s/,//g;
    $cumulative_n_proteins =~ s/,//g;

		if ($idx ==0 ){
				push @cumpepx, 0;
			}else{
				push @cumpepx, $pre_cum_n_good_spectra;
        push @cumpepx, $pre_cum_n_good_spectra;
			}
			push @cumpepx, $pre_cum_n_good_spectra+$n_good_spectra;
			push @cumpepy, ($cumulative_n_peptides,$cumulative_n_peptides,0);
			push @idvpepy, ($n_distinct_peptides,$n_distinct_peptides,0);
			push @cumproty, ($cumulative_n_proteins,$cumulative_n_proteins,0);
			push @idvproty, ($n_canonical_proteins,$n_canonical_proteins,0);
			$pre_cum_n_good_spectra += $n_good_spectra;
			$idx++;
  }
  my $cumpepx_str = join(",", @cumpepx);
  my $cumpepy_str = join(",", @cumpepy);  
  my $idvpepy_str = join(",", @idvpepy);
  my $cumproty_str = join(",", @cumproty);
  my $idvproty_str = join(",", @idvproty);
  my $sample_label_str = join("','", @sample_label);
  my $total_spec = $cumpepx[$#cumpepx];
  my $width =  $total_spec < 5000000 ? 'width: 1000,' : '';
  my $plot_js = qq~
			l=['$sample_label_str']
			var cum = {
				x: [$cumpepx_str],
				y: [$cumpepy_str],
				fill: 'tonexty',
        fillcolor:'#1f77b4',
        marker: {
                  color: 'rgba(255,0,255,0)',
                }, 
				name:'cumulative_n_peptides',
				hovertext:l,
				haveron:"fills",
        hoverinfo:"x+y+text",
			};
			var idv = {
				x: [$cumpepx_str],
				y: [$idvpepy_str],
				fill: 'tonexty',
        fillcolor:'#ff7f0e',
				marker: {
					color: 'rgba(255,0,255,0)',
				},
				name:'n_distinct_peptides',
				hovertext:l,
        hoverinfo:"x+y+text",
				haveron:"fills"
			};

			var layout = {
				$width
				height: 800,
				font: {
					size: 18
				},
        legend:{x: 0.029,y: 1.1,font: { size: 12}},
				xaxis:{title:'Cumulative Number of MS/MS Spectra Identified'},
				yaxis:{title:'Number of Distinct Peptides'}
			};
			var data = [idv,cum];
			var config = {
				toImageButtonOptions: {
					format: 'png', // one of png, svg, jpeg, webp
					filename: 'image',
					height: 800,
					width: 1000,
					scale: 6,
				}
			};
			Plotly.newPlot('plot_div', data,layout,config);

      var cum2 = {
        x: [$cumpepx_str],
        y: [$cumproty_str],
        fill: 'tonexty',
        fillcolor:'#003f72',
        marker: {
                  color: 'rgba(255,0,255,0)',
                },
        name:'cumulative_n_proteins',
        hovertext:l,
        hoverinfo:"x+y+text",
        haveron:"fills"
      };
      var idv2 = {
        x: [$cumpepx_str],
        y: [$idvproty_str],
        fill: 'tonexty',
        fillcolor:'#b00',
        marker: {
                  color: 'rgba(255,0,255,0)',
                },
        name:'n_canonical_proteins',
        hovertext:l,
        hoverinfo:"x+y+text",
        haveron:"fills"
      };
      var layout = {
         $width
        height: 800,
        font: {
          size: 18
        },
        legend:{x: 0.029,y: 1.1,font: { size: 12}},
        xaxis:{title:'Cumulative Number of MS/MS Spectra Identified'},
        yaxis:{title:'Number of Distinct Proteins'}
      };
      var data2 = [idv2,cum2];
      Plotly.newPlot('plot_div2', data2,layout,config);
  ~;
  #print "$plot_js<BR>";

  my $chart = qq~
     <!-- Latest compiled and minified plotly.js JavaScript -->
     <script type="text/javascript" src="https://cdn.plot.ly/plotly-latest.min.js"></script>
     <p class="plot_caption"><b>Plot below shows the number of peptides contributed by each experiment</b>, and the cumulative number of distinct peptides for the build as of that experiment.</p> 
     <div style="width: 80vw"><div id="plot_div"></div><br></div>
     <p class="plot_caption"><b>Plot below shows cumulative number of canonical proteins contributed by each experiment.</b><br>
     Height of red bar is the number of proteins identified in experiment; 
     height of blue bar is the cumulative number of proteins;<br> 
     width of the bar (x-axis) shows the number of spectra identified (PSMs), above the threshold, for each experiment.</p>
     <div style="width: 80vw"><div id="plot_div2" ></div></div><br><br>
	<script type="text/javascript" charset="utf-8">
      $plot_js
	</script>
  ~;
  return $chart;
  
}
sub tableHeatMap{
  my $self = shift;
  my %args = @_;
  my $table_id = $args{table_id};
  my $col = $args{column} || 2;
  my $total = $args{total} || 0;
  my $str = qq~
     <script type="text/javascript" src="$CGI_BASE_DIR/../usr/javascript/jquery/jquery.min.js"></script>
     <script type="text/javascript">
      jQuery.noConflict();
      jQuery(document).ready(function(){
	  // Function to get the Max value in Array
	  Array.sum = function( array ){
            var sum =0;
	    for (let i in array){
              sum = sum+ array[i];
            }
            return sum;
	  };

	  // get all values
	  var counts= jQuery('#$table_id  td:nth-child($col)').map(function() {
	    if (jQuery(this).text() ){
	      return parseInt(jQuery(this).text());
	    }else{
	      return 0;
	    }
	}).get();
	// return max value
	  var sum = Array.sum(counts);
        if ($total > 0){
           sum = $total;
        }
	// add classes to cells based on nearest 10 value
	jQuery('#$table_id  td:nth-child($col)').each(function(){
		var val = parseInt(jQuery(this).text());
		var pctval = parseInt((Math.round((val/sum)*100)).toFixed(0));
                var pctval2 = 100 - pctval ;
		clr = "linear-gradient(to right, #f46d69 " + pctval + "%, #ffffff00 " + pctval + "% " + pctval2 + "%)"; 
		jQuery(this).css("background-image", clr);
	     });
	});
	</script>
  ~;
  return $str;
}
##################################################################################
sub plotly_pie {

  my $self = shift;
  my %args = @_;
  my $data = $args{data};
  my $names = $args{names};
  my $divname = $args{divName} || 'pie_div';
  
  my $values = join(",", @$data);
  my $labels = join('","', @$names);
    
  ###<script src="https://cdn.plot.ly/plotly-latest.min.js"></script> 
	my $plot_js = qq~
   <script type="text/javascript" charset="utf-8">
   var data=[{
      values: [$values],
      labels:["$labels"],
      type: 'pie',
      //textinfo: "label+percent",

    }]
		var pie_layout = {
			height: 400,
			width: 800
		};
    Plotly.newPlot("$divname", data, pie_layout);
	 </script>
  ~;
  return $plot_js;
}

##################################################################################
sub plotly_barchart {
  my $self = shift;
  my %args = @_;
  my $all_data = $args{data};
  my $names = $args{names};
  my $divname = $args{divName} || '';
  my $xtitle = $args{xtitle} || '';
  my $ytitle = $args{ytitle} || '';
  my $title = $args{title} || '';
  my $dtick = $args{dtick} || 5;
  my $colors = $args{colors} ;
  my $xtickangle = $args{xtickangle} || '';
  my $layoutmargin = $args{layoutmargin} || '';
  my $yticktype = $args{yticktype} || '';
  my $barlabel = $args{barlabel} || 0;

  my $plot_js = '';
  my $counter = 0;
  my @ts = ();
  my $layout_title = "title:'<b>$title</b>'";
  my $n_series = scalar @$names;
  my $bargap = ",'bargap':0.7" if ($n_series == 1);

  foreach my $data (@$all_data){
		my @category = () ;
		my @cnt = () ;
		foreach my $row (@$data){
			push @category, $row->[0];
			push @cnt , $row->[1];
		}
		my $category_str = join("','", @category);
		my $cnt_str = join(",", @cnt);
    my $marker = '';
    my $bartext = '';

    if ($colors && $colors->[$counter]){
      $marker = "marker: {color: '$colors->[$counter]";
    }
    if ($barlabel){
       $bartext = "text:y$counter.map(String),textposition:'outside',cliponaxis:false,";
    }
		$plot_js .= qq~
        y$counter = [$cnt_str];
				var t$counter = {
					x: ['$category_str'],
					y: y$counter,
					name: '$names->[$counter]', 
					type: 'bar',	
          hoverinfo: 'text',
          $bartext
				  $marker	
				};
			~;
     push @ts, "t$counter";
     $counter++;
   }
   my $ts_str = join(",", @ts);
   $plot_js .= qq~
				var data = [$ts_str];
				var config = {
					toImageButtonOptions: {
						format: 'png', // one of png, svg, jpeg, webp
						filename: 'image',
						height: 800,
						width: 1000,
						scale: 6,
					}
				};
				var layout = {
					barmode: 'group',
          font: {size: 12},
          //legend:{x: 0.02,y: 1.25 ,font: { size: 12}},
          legend: {x:1.1,  xanchor: 'right',y:1, font: { size: 12}},
          xaxis: {type: 'category',
                  dtick:$dtick,
                  title: '$xtitle', 
                  ticklabeloverflow:'allow',
                  $xtickangle
                   },
          yaxis: {title: '$ytitle',
                  $yticktype
                 },
          $layout_title,
          hovermode: 'closest',
          margin:{$layoutmargin}
          $bargap
        };
        // Use 'closest' to show only the hover label without marker
        Plotly.newPlot('$divname', data, layout,config);
   ~;
  my $chart = qq~
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script> 
		<script type="text/javascript" charset="utf-8">
      $plot_js
		</script>

  ~;
  return $chart;
}

sub get_dataset_url{
  my $self = shift;
  my $repository_ids = shift;
  my $url ='';
  my @ids = split(/[,;]/, $repository_ids);
	my $annotation_url = '';
  unless (%dataset_annotation){
		my $url = "http://proteomecentral.proteomexchange.org/api/autocomplete/v0.1/datasets";
		my $content  = get($url);
		$content =~ s/[\[\]"\s+\n\r]//g;
		foreach my $pxd (split(",", $content)){
			next if ($pxd eq '');
			$dataset_annotation{$pxd} =1;
		}
  }
  $url='';

  foreach my $id(@ids){
    $id =~ s/\s+//g;
    my $link;
    if ($id =~ /PXD/){
      $link = "http://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=";
    }elsif($id =~ /PASS/){
      $link= "http://www.peptideatlas.org/PASS/";
    }elsif($id =~ /^S\d+/){
      $link= "https://cptac-data-portal.georgetown.edu/study-summary/";
    }
    if ($link){
      $url .= "<a href='$link$id' target='_blank'>$id</a>";
      if (defined $dataset_annotation{$id}){
        $annotation_url = "<a href='http://www.peptideatlas.org/datasets/annotation.php?dataset_id=$id' target='_blank'>annot</a>";
      }
      $url .= ",";
    }else{
      $url .= "$id,";
    }
  }
  $url =~ s/,$//;
  return ($url,$annotation_url);
}

sub create_table {
  my $self = shift;
  my %args = @_;
  my $data = $args{data} || die "need data\n";
  my $column_names = $args{column_names} || die "need column_names\n";
  my $table_name =  $args{table_name} || die "need table name\n";
  my $table_id =  $args{table_id} || die "need div name\n";
  my $plot_html = $args{plot_html} || '';
  my $align = $args{align} || ();
  my $sortable = $args{sortable} || 0;
  my $table_width = $args{width} || 800;
  my $rows_to_show = $args{rows_to_show} || 25;
  my @noWrapCol = ();
  my $nowrap = $args{nowrap} || \@noWrapCol; 
  my $download_table = $args{download_table} || 0;
  my $default_sort_col = $args{default}  || '';
  if ($sortable){
    my @sort_headings = ();
    foreach my $col (@$column_names){
       push @sort_headings, $col, '';
    }
    
    unshift @$data, $self->make_sort_headings( headings => \@sort_headings, default => $default_sort_col );
  }else{
    unshift @$data, $column_names;
  }
 
  my $table = $self->encodeSectionTable(unified_widgets => 1,
                              rows => $data,
                              header => 1,
															header_sticky => 1,
                              bkg_interval => 3,
                              set_download => $download_table,
                              nowrap => $nowrap,
                              table_id => $table_id,
                              align => [@$align],
                              width => $table_width ,
                              rows_to_show => $rows_to_show, 
                              sortable => $sortable );


  my $heading_info = $self->get_table_help(column_titles_ref=> $column_names);
  my $html = $sbeams->make_toggle_section(
      neutraltext =>$table_name,
      sticky => 1,
      barlink => 1,
      visible => 1,
      name => $table_id."_div",
      content =>"<TABLE><TR><TD COLSPAN='5'>$plot_html</td></tr><TR><TD COLSPAN='5'>$heading_info</TD></TR>$table</TABLE>",
  );

  return $html;
}

sub display_spectra_ptm_table {
  my $self = shift;
  my %args = @_;
  my $ptm_identification_ref = $args{ptm_identification_ref};
  my $column_titles_ref = $args{column_titles_ref};
  my $colnameidx_ref = $args{colnameidx_ref};
  my $hidden_cols_ref= $args{hidden_cols_ref};
  my $resultset_ref = $args{resultset_ref};
  my $ptm_score_summary_ref = $args{ptm_score_summary_ref};
  my $drawVisualization = '';

  my $spectraHTML = qq~
    <div class="tab">
      <style>
      /* Style the tab */
      .tab {
        overflow: hidden;
      }
      /* Style the buttons inside the tab */
      .tab button {
        padding: 14px 16px;
        border: none;
        outline: none;
        font-size: 14px;
        font-weight:bold;
        float:left;
        color:#555555;
        text-decoration:none;
        background:#f3f1e4;
        border-right:1.25px solid #AAAAAA;
      }

      /* Change background color of buttons on hover */
      .tab button:hover {
        background-color:#bb0000;
            color:#ffffff;
      }

      /* Create an active/current tablink class */
      .tab button.active {
         background:#ffffff;
          color:#333333;
      }

      /* Style the tab content */
      .tabcontent {
        border-top: none;
      }
		 </style>
		 ~;
	#$spectraHTML = $sbeams->getPopupDHTML();
	my $counter = 1;
	foreach my $ptm_type (keys %$ptm_identification_ref){
		$spectraHTML .="<button type='button' class='ptmtablinks' onclick='openPTM(event,\"$ptm_type\")'>$ptm_type</button>\n";
		$counter++;
	}
  $spectraHTML .="</div>";

  $counter = 1;
  foreach my $ptm_type(keys %$ptm_identification_ref){
    my @data =();
    foreach my $s ( sort {$b <=> $a} keys %{$ptm_identification_ref->{$ptm_type}}){
      push @data, @{$ptm_identification_ref->{$ptm_type}{$s}};
    }

		$resultset_ref->{data_ref} = \@data;
    $hidden_cols_ref->{'ptm_lability'} = 0;
    $hidden_cols_ref->{'ptm_sequence'} = 0;
    #$hidden_cols_ref->{'ptm_norm_info_gain'} = 0;
		my $spectra_display = $self->get_individual_spectra_display(
			column_titles_ref=>$column_titles_ref,
			colnameidx_ref => $colnameidx_ref,
			resultset_ref=> $resultset_ref,
			hidden_cols_ref=>$hidden_cols_ref,
      table_id => 'ptm_spectra' );
		$spectra_display =~ s/<table/<table name="ptm_spectra_$counter" class="ptm_spectra"/i;

		my $GV = SBEAMS::Connection::GoogleVisualization->new();
		my $ptm_summary_chart = $GV->drawPTMHisChart(data=> $ptm_score_summary_ref->{$ptm_type},ptm_type=>$ptm_type );
    my $func_name = "drawVisualization_$ptm_type";
    $func_name =~ s/[:\.\-\+\(\)]/_/g;
		if ($counter == 1){
			$spectraHTML .= "<div  id ='$ptm_type' style='display:block' class='tabcontent'>\n";
		}else{
			$spectraHTML .= "<div  id ='$ptm_type' style='display:none' class='tabcontent'>\n";
      $drawVisualization .= "google.setOnLoadCallback($func_name)\n";
		}
		$spectraHTML .= qq~<div>$ptm_summary_chart</div>\n~;
		#$spectraHTML .= "$ptm_prob_color\n$spectra_display";
		$spectraHTML .= "$spectra_display\n</div>";
		$counter++;
  }
		#### Display table
	print scalar $sbeams->make_toggle_section( neutraltext => "PTM Spectra",
							 sticky => 1,
							 name => 'ptm_spectra_div',
							 barlink => 1,
							 visible => 1,
							 content => $spectraHTML );

		print qq~
		<script LANGUAGE="JavaScript" TYPE="text/javascript">
     function openPTM(evt, Name) {
        $drawVisualization;     
        var i, tabcontent, ptmtablinks;
        tabcontent = document.getElementsByClassName("tabcontent");
        for (i = 0; i < tabcontent.length; i++) {
          tabcontent[i].style.display = "none";
        }
        ptmtablinks = document.getElementsByClassName("ptmtablinks");
        for (i = 0; i < ptmtablinks.length; i++) {
          ptmtablinks[i].className = ptmtablinks[i].className.replace(" active", "");
        }
        document.getElementById(Name).style.display = "block";
        evt.currentTarget.className += " active";
      }

			var indiv_spec = document.getElementsByClassName("ptm_spectra");
      for (var k=0; k<indiv_spec.length; k++) {
				var rows = indiv_spec[k].getElementsByTagName('TR');
				var j=3;
      }
		</script>
  ~;

}
sub get_pfam_domain_display {
  my $self = shift;
  my %args = @_;
  my $data = $args{data}; ## json format
  my $list = $args{list};
  my $html = qq~
    <script type="text/javascript" src="$HTML_BASE_DIR/usr/javascript/pfam/prototip/prototype.js"></script>
    <script type="text/javascript" src="$HTML_BASE_DIR/usr/javascript/pfam/prototip/prototip.js"></script>
    <script type="text/javascript" src="$HTML_BASE_DIR/usr/javascript/pfam/prototip//styles.js"></script>
    <link rel="stylesheet" href="$HTML_BASE_DIR/usr/javascript/pfam/prototip/prototip.css" type="text/css" />

    <!-- stylesheets. We only really need the rules that are specific to the tooltips -->
    <link rel="stylesheet" href="$HTML_BASE_DIR/usr/javascript/pfam/css/pfam.css" type="text/css" />

    <script type="text/javascript" src="$HTML_BASE_DIR/usr/javascript/pfam/js/main.js"></script>

    <script type="text/javascript" src="$HTML_BASE_DIR/usr/javascript/pfam/js/excanvas.js"></script>
    <script type="text/javascript" src="$HTML_BASE_DIR/usr/javascript/pfam/js/graphic_generator.js"></script>
    <script type="text/javascript" src="$HTML_BASE_DIR/usr/javascript/pfam/js/domain_graphics.js"></script>

  ~;

  my $counter = 1;
  $html .= '<div id="graphic" class="panel">';
  $html .='<h2>Domain graphic</h2>';

  foreach my $prot (@$list){
    $html .=qq~
    <div style="display:flex;width:100%;text-align: center;">
    <div style="float:left;width:15%;text-align:right;flex-grow:1;">$prot</div>
    ~;
    if(defined $data->{$prot}){
		  $html .= qq~<div style="float:right;width:85%;flex-grow:1;text-align:left;" id="dg$counter"><span id="none"></span></div>~;
      $counter++;
    }else{
      $html .= qq~<div style="float:right;width:85%;flex-grow:1;text-align:left;"></div>~;
    }
    $html .="</div>";
  }
  $html .=qq~
	</div>
  <br>
  <div class="cleaner"><!-- empty --></div>
	<div id="form">
	<div id="sequence" class="panel">
  ~;
  
  use JSON;
  my $json = JSON->new->allow_nonref;;
  $counter=1;
  foreach my $prot (@$list){
		$json = $json->pretty([1]);
		my $json_text =   $json->encode($data->{$prot});
		$html .= "<textarea cols='80' rows='40' id='pfam_seq$counter' hidden>$json_text</textarea>\n";
    $counter++;
  }
  $html .= qq~
		</div>
  	</div>

	<div> 
		<div class="row">
		<span class="label">X-scale:</span><input id="xscale" value="2.0"></input>
		</div>
		<div class="row">
		<span class="label">Y-scale:</span><input id="yscale" value="1.0"></input>
		</div>
		<button id="submit">Update graphic</button>
		<div id="errors" style="display: none"></div>
		<div id="error"></div>
	</div>

    <script type="text/javascript">
  var generator;
  document.observe( "dom:loaded", function() {
    generator = new GraphicGenerator();
  } );

    </script>
  ~;
   
  return $html;
   


}
sub draw_tree{

	my $self = shift;
	my $newick = shift;
	my $html = qq~ 
    <script src="https://code.jquery.com/jquery.js"></script>
    <script src="$HTML_BASE_DIR/usr/javascript/tree/underscore-min.js"></script>
    <link rel="stylesheet" href="https://netdna.bootstrapcdn.com/font-awesome/4.0.3/css/font-awesome.css"/>
    <script src="$HTML_BASE_DIR/usr/javascript/tree/bootstrap.min.js" ></script>
    <script src="$HTML_BASE_DIR/usr/javascript/tree/d3.v5.js"></script>
    <script src="$HTML_BASE_DIR/usr/javascript/tree/phylotree.js"></script>
    <link href="$HTML_BASE_DIR/usr/javascript/tree/phylotree.css" rel="stylesheet" />

    <style>

      \@media (max-width: 1075px) {
        .container {
          padding-top: 50px;
        }
      }

      .btn {
        height: 30px;
      }

      #toolbar {
        display: flex;
        justify-content: space-between;
        width: 550px;
      }

      #controls {
        position: fixed;
        left: 5px;
      }
    </style>

    <div class="container" id="main_display">
      <div class="row">
        <div class="col-md-8">
          <div class="btn-toolbar" role="toolbar" id="toolbar">
            <div class="btn-group">
              <button
                type="button"
                class="btn btn-light btn-sm"
                data-direction="vertical"
                data-amount="1"
                title="Expand vertical spacing"
              >
                <i class="fa fa-arrows-v"></i>
              </button>
              <button
                type="button"
                class="btn btn-light btn-sm"
                data-direction="vertical"
                data-amount="-1"
                title="Compress vertical spacing"
              >
                <i class="fa  fa-compress fa-rotate-135"></i>
              </button>
              <button
                type="button"
                class="btn btn-light btn-sm"
                data-direction="horizontal"
                data-amount="1"
                title="Expand horizonal spacing"
              >
                <i class="fa fa-arrows-h"></i>
              </button>
              <button
                type="button"
                class="btn btn-light btn-sm"
                data-direction="horizontal"
                data-amount="-1"
                title="Compress horizonal spacing"
              >
                <i class="fa  fa-compress fa-rotate-45"></i>
              </button>
              <button
                type="button"
                class="btn btn-light btn-sm"
                id="sort_ascending"
                title="Sort deepest clades to the bototm"
              >
                <i class="fa fa-sort-amount-asc"></i>
              </button>
              <button
                type="button"
                class="btn btn-light btn-sm"
                id="sort_descending"
                title="Sort deepsest clades to the top"
              >
                <i class="fa fa-sort-amount-desc"></i>
              </button>
              <button
                type="button"
                class="btn btn-light btn-sm"
                id="sort_original"
                title="Restore original order"
              >
                <i class="fa fa-sort"></i>
              </button>
              <button
                type="button"
                class="btn btn-light btn-sm"
                id="save_image"
                title="Save image"
              >
                <i class="fa fa-picture-o"></i>
              </button>
            </div>
            <div class="btn-group" role="group">
              <button
                class="btn btn-light btn-sm active phylotree-layout-mode"
                data-mode="linear"
              >
                Linear
              </button>
              <button
                class="btn btn-light btn-sm phylotree-layout-mode"
                data-mode="radial"
              >
                Radial
              </button>
            </div>
            <div class="btn-group" role="group">
              <button
                class="btn btn-light btn-sm active phylotree-align-toggler"
                data-align="left"
              >
                <i class="fa fa-align-left"></i>
              </button>
              <button
                class="btn btn-light btn-sm phylotree-align-toggler"
                data-align="right"
              >
                <i class="fa fa-align-right"></i>
              </button>
            </div>
          </div>
        </div>

      </div>

      <div class="row">
        <div class="col-md-12">
          <div id="tree_container" class="tree-widget"></div>
        </div>
      </div>
    </div>
    <BR>
        <script src="../../usr/javascript/tree/plot_phylotree.js"></script>
    <script>
      var newick_string = '$newick';
    </script>
  ~;
  return $html;
}

sub draw_venn {
  my $self = shift;
  my %args =@_;
  my @data = @{$args{data}};
  my $data_label = $args{data_label} || ();
  my $title = $args{title} || '';
  my $divid = $args{id} || 'venn';
  my $proportion = $args{proportion} || 'yes';

  return '' if (! @data || @data == 1 );
  my %values = ();
  my %summary = ();
  my %unique=();
  my $n_set = scalar @data ;
  my $i=1;
  my @n = (0..$n_set-1);

  sub combine;
  while ($i <= $n_set){
    foreach (combine [@n], $i){
      $summary{join(",", @$_)} = 0
    } 
    $i++;
  }

  for (my $i=0; $i<$n_set; $i++){
    my $n=scalar @{$data[$i]};
    $summary{$i} = $n;
    foreach my $val (@{$data[$i]}){
			$values{$val}{$i}=1;
    }
  }

  foreach my $val (keys %values){
    my $set  = join(',', sort {$a <=> $b} keys %{$values{$val}});
    if ($set =~ /,/){
      $summary{$set}++;
    }else{
      $unique{$set}++;
    }
  }

  my $html = qq~
		<script src="../../usr/javascript/d3_venn/d3.v4.min.js"></script>
		<script src="../../usr/javascript/d3_venn/venn.js"></script>
    <script> 
    var sets=[~;

  my $sep = ''; 
  foreach my $s (sort {$a <=> $b} keys %summary){
    my $text = '';
    if ($s !~ /,/){
      $text = "$data_label->[$s] ". ($unique{$s} || 0);
    }else{
      $text ="$summary{$s}";
    }
    if( $proportion eq 'no'){
      my $size = 100;
      if (length($s) > 1){
        my @m = $s =~ /(\d+)/g;
        my $n = scalar @m;
        $size = ($n-1)/($n_set-1) * 50;
      }
      $html .= qq~$sep\n{"sets": [$s], "label": "$text", "size": $size}~;
    }else{
      $html .= qq~$sep\n{"sets": [$s], "label": "$text", "size":$summary{$s}}~;
    }
    $sep = ',';
  }
  $html .= "];\n";
  $html .=qq~
      var chart = venn.VennDiagram()
                             .width(300)
                             .height(300);
			d3.select("#$divid").datum(sets).call(chart);
      var areas = d3.selectAll("#$divid g");
      areas.select("text")
           .style("fill-opacity", 1)
					.style("font-weight", "bold");

   </script>
  ~;

}

sub combine {

  my ($list, $n) = @_;
  die "Insufficient list members" if $n > @$list;
  return map [$_], @$list if $n <= 1;
  my @comb;
  for (my $i = 0; $i+$n <= @$list; ++$i) {
    my $val  = $list->[$i];
    my @rest = @$list[$i+1..$#$list];
    push @comb, [$val, @$_] for combine \@rest, $n-1;
  }
  return @comb;
}

1;

__END__
###############################################################################
###############################################################################
###############################################################################

=head1 NAME

SBEAMS::WebInterface::HTMLPrinter - Perl extension for common HTML printing methods

=head1 SYNOPSIS

  Used as part of this system

    use SBEAMS::WebInterface;
    $adb = new SBEAMS::WebInterface;

    $adb->printPageHeader();

    $adb->printPageFooter();

    $adb->getGoBackButton();

    $adb->getGoBackButton();

=head1 DESCRIPTION

    This module is inherited by the SBEAMS::WebInterface module,
    although it can be used on its own.  Its main function 
    is to encapsulate common HTML printing routines used by
    this application.

=head1 METHODS

=item B<printPageHeader()>

    Prints the common HTML header used by all HTML pages generated 
    by theis application

=item B<printPageFooter()>

    Prints the common HTML footer used by all HTML pages generated 
    by this application

=item B<getGoBackButton()>

    Returns a form button, coded with javascript, so that when it 
    is clicked the user is returned to the previous page in the 
    browser history.

=head1 AUTHOR

Eric Deutsch <edeutsch@systemsbiology.org>

=head1 SEE ALSO

perl(1).

=cut


