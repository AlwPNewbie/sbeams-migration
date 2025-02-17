#!/usr/local/bin/perl -w

use strict;
use Storable qw(nstore retrieve);
use Test::More tests => 6;
use Test::Harness;
use POSIX;

my $aname = POSIX::tmpnam();
ok( store_array(), 'Store array');
ok( retrieve_array(), 'Retrieve array');

my $hname = POSIX::tmpnam( );
ok( store_hash(), 'Store hash');
ok( retrieve_hash(), 'Retrieve hash');

my $rname = POSIX::tmpnam( );
ok( retrieve_rs(), 'Retrieve resultset');
SKIP: {
        skip "FIXME!", 1;
ok( retrieve_DATA_rs(), 'Retrieve DATA resultset');
      }


sub store_array {
  my @array = ( 'PeptideAtlas_atlas_build_id', 123 );
  nstore(\@array, $aname ) || return 0;
}

sub store_hash {
  my %hash = ( PeptideAtlas_atlas_build_id => 123 );
  nstore(\%hash, $hname ) || return 0;
}

sub retrieve_array {
  my $arrayref = retrieve( $aname ) || return 0;
  my $ok = $arrayref->[0];
  return $ok
}

sub retrieve_hash {
  my $hashref = retrieve( $hname ) || return 0;
  my $ok = 0;
  for my $k ( keys( %$hashref ) ) {
    $ok = $k;
  }
  unless ( -e 'test.sto' ) {
    print STDERR "$hname doesn't exist!\n";
    `cp $hname test.sto`;
  }
  return $ok
}

sub retrieve_DATA_rs {
  open( RS, ">$rname" ) || die "unable to open $rname";
  print STDERR "$rname\n";

  while( my $rs = <DATA> ) {
    print RS $rs;
	}
  close RS;
  my $ok = 0;

  my $hashref = eval {retrieve( $rname )} || return 0;
  for my $k ( keys( %$hashref ) ) {
    $ok = $k;
  }
  return $ok
}

sub retrieve_rs {
  my $hashref = retrieve( 'test.sto' ) || return 0;
  my $ok = 0;
  for my $k ( keys( %$hashref ) ) {
    $ok = $k;
  }
  return $ok
}

__DATA__
pst01234      
int
12
12
12
12   types_list_ref�   row_counter   

project_id
project_tag
name
username
description   column_list_ref�   row_pointer   �   project_tag�   name�   description�   username�
   project_id   column_hash_ref�	   page_size   �����   precisions_list_ref       data_ref
