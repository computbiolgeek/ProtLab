#!/usr/bin/perl

use strict;
use warnings;

use Getopt::Long qw(GetOptions);
use Carp qw(carp croak);
use Data::Dumper;

my $db;
my $pdb_list;
my $output;

# set up commandline options
GetOptions(
	"pdb_list=s" => \$pdb_list,
	"output=s"   => \$output,
	"db=s"       => \$db,
) or die usage();

# check whether required commandline paramemters are provided
unless ( defined $pdb_list && defined $output && $db ) {
	die usage();
}

unless ( -d $db ) {
	carp "$db not exists or is not a directory, create $db ...";
	mkdir($db);
}

# the following is the main logic of the task
open my $fh, "<", $pdb_list or die "Could not open $pdb_list for reading: $!";
chomp( my @pdb_ids = <$fh> );
close($fh);

# open file for writing
open my $ofh, ">", $output or die "Could not open $output for writing: $!";
foreach my $pdb_id (@pdb_ids) {
	get_pdb_from_opm( $pdb_id, $db );
	my $pdb_file       = $db . "/" . $pdb_id . ".pdb";
	my %chain2cacoords = get_chain2cacoords_hash($pdb_file);
	my @tm_chain_ids   = get_tm_chains( \%chain2cacoords );
	unless (@tm_chain_ids) {
		carp "No transmembrane subunits found in $pdb_id, please verify";
	}
	foreach my $id (@tm_chain_ids) {
		print $ofh $pdb_id . $id . "\n";
	}
}
close($ofh);

######################INTERFACE SUBROUTINE##########################
# Comments     : For each chain in the PDB file, store the coordinates 
#                of all CAs in a hash of arrays.
# Usage        : get_chain2cacoords_hash($pdb_file)
# Parameter(s) : the pdb file to be processed
# Returns      : a hash of chain ID to array of CA coordinates
####################################################################
sub get_chain2cacoords_hash {

	# pdb file name
	my ($pdb_file) = shift @_;

	# parse the pdb file for relevant information
	print "Parsing $pdb_file for membrane half thickness and CA coordinates\n";

	# extract half membrane thickness
	# this information is on the first line of the opm pdb file
	open my $pdb_fh, "<", $pdb_file
	  or die "Could not open $pdb_file for reading: $!";
	my $first_line = <$pdb_fh>;
	chomp($first_line);
	my @fl_tokens = split( /\s+/, $first_line );
	my $half_thickness = $fl_tokens[-1];

	# read the rest records into an array, one record at a time
	chomp( my @pdb_records = <$pdb_fh> );
	my @ca_records =
	  grep { $_ =~ /^ATOM/ && substr( $_, 12, 4 ) eq " CA " } @pdb_records;
	close($pdb_fh);

	# hash: chain id to CA coordinates of the corresponding chain
	my %chain2cacoords = ( HALF_THICKNESS => $half_thickness, );

	# populate the hash
	foreach my $record (@ca_records) {
		my $chainid = substr( $record, 21, 1 );

		# x, y, z coordinates are currently stored as a string
		# which will be split in get_tm_chains
		my $ca_coordinates = substr( $record, 30, 24 );
		push( @{ $chain2cacoords{$chainid} }, $ca_coordinates );
	}

	# return the hash
	return %chain2cacoords;
}

######################INTERFACE SUBROUTINE##########################
# Comments     : Given a hash of chain IDs to CA coordinations, returns
#                the IDs of chains that are transmembrane.
# Usage        : get_tm_chains($hashref)
# Parameter(s) : a reference to a hash of chain IDs to CA coordinations
# Returns      : the IDs of chains that are transmembrane
####################################################################
sub get_tm_chains {
	my ($chain2cacoords) = shift @_;
	my @tm_chain_ids;
	my $half_thickness = $chain2cacoords->{HALF_THICKNESS};
	delete $chain2cacoords->{HALF_THICKNESS};
	while ( my ( $key, $ca_coords ) = each %{$chain2cacoords} ) {
		if ( is_tm_chain( $half_thickness, $ca_coords ) ) {
			push( @tm_chain_ids, $key );
		}
	}
	return @tm_chain_ids;
}

######################INTERFACE SUBROUTINE##########################
# Comments     : Given the half-thickness of the membrane and all
#                the CA coordiates of a chain, determine whether this
#                chain is transmembrane.
# Usage        : is_tm_chain($half_thickness, $ca_coords)
# Parameter(s) : 
#                $half_thickness: the half-thickness of the membrane
#                $all_ca_coords : reference to an array of CA coords
# Returns      : boolean value indicating whether the chain is TM
####################################################################
sub is_tm_chain {
	my ( $half_thickness, $all_ca_coords ) = @_;
	my $find_residue_above_membrane = 0;
	my $find_residue_below_membrane = 0;

	# iterate over all CA atoms to find one above the membrane
	# and one below the membrane
	foreach my $coord_str ( @{$all_ca_coords} ) {

		# get z coordinate
		my $z = substr($coord_str, 16, 8);
		if ( !$find_residue_above_membrane && $z > $half_thickness ) {
			$find_residue_above_membrane = 1;
		}
		elsif ( !$find_residue_below_membrane && $z < -$half_thickness )
		{
			$find_residue_below_membrane = 1;
		}
		if ( $find_residue_above_membrane && $find_residue_below_membrane ) {
			return 1;
		}
	}
	return 0;
}

######################INTERFACE SUBROUTINE##########################
# Comments     : Retrieve the requested PDB file from OPM database.
# Usage        : get_pdb_from_opm($pdb_id, $db)
# Parameter(s) :
#                $pdb_id: the requested PDB ID
#                $db    : path to the directory that stores the 
#                         requested PDB file
# Returns      : 
####################################################################
sub get_pdb_from_opm {
	my ( $pdb_id, $db ) = @_;
	my $pdb_filename = $db . "/" . $pdb_id . ".pdb";
	if ( !-f $pdb_filename ) {
		print "Downloading $pdb_id.pdb from the OPM database\n";
		system( "curl -m 60 -s http://opm.phar.umich.edu/pdb/"
			  . $pdb_id
			  . ".pdb -o "
			  . $pdb_filename );
	}
	print "$pdb_filename exits, skip\n";
}

######################INTERFACE SUBROUTINE##########################
# Comments     : Print help infomation.
# Usage        : usage()
# Parameter(s) : 
# Returns      : 
####################################################################
sub usage {
	print "Usage: $0 --pdb_list <a file containing a list of pdb ids> --output <output file name>\n";
}

=head1 NAME

extract_tm_chains.pl - Given a list of PDB IDs, this scripts extract transmembrane chains.

=head1 VERSION

This documentation refers to extract_tm_chains.pl version 0.0.1.

=head1 DESCRIPTION

This script takes a list of pdb ids and extracts chains in each pdb that are transmembrane.

=head1 USAGE

extract_tm_chains.pl --pdb_list pdb_list.txt --db ~/database --output tm_chains.txt

=head1 AUTHOR

Bian Li, bian.li@vanderbilt.edu

=head1 COPYRIGHT

Copyright (c) 2016 Bian Li. All rights reserved.

=head1 LICENSE

This application is free software; you can redistribute it and/or modify it under
the same terms as Perl itself.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.