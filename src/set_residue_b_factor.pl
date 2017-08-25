#!/usr/bin/perl

####################################################
# Sets the b-factor column of a pdb file
# to user-provided b-factors.
# @author Bian Li
# Please send bug information to Bian Li
# via bian.li@vanderbilt.edu
####################################################

use strict;
use warnings;

unless(@ARGV) {
	Usage();
	exit;
}

# global variables
my $b_factor_file="";
my $input_pdb="";
my $output_pdb="";

# parse commandline arguments
ParseInput();

# read the b-factor file and the input pdb file
open(B_FACTOR, "< $b_factor_file") or die "Could not open $b_factor_file for reading $!";
my @b_factors = ();
while(my $line = <B_FACTOR>) {
	chomp($line);
	# skip uninterested lines
	if($line =~ /^\s*#|^\s*$/) {
		next;
	}
	my($index, $value) = ($line =~ /(\S+)\s*(\S+)/);
	$b_factors[$index] = $value; 
}
close(B_FACTOR);

open(INPUT_PDB, "< $input_pdb") or die "Could not open $input_pdb for reading $!";
my @pdb_lines = <INPUT_PDB>;
close(INPUT_PDB);

# open output pdb file for writing
open(OUTPUT_PDB, "> $output_pdb") or die "Could not open $output_pdb for writing";

# change the b-factor column
foreach my $line (@pdb_lines) {
	chomp($line);
	if($line =~ /^ATOM/) {
		my $residue_id = substr($line, 22, 4);
		$residue_id =~ s/^\s*//;
		my $b_factor = 0;
		if(exists $b_factors[$residue_id]) {
			$b_factor = $b_factors[$residue_id];
		}
		$line = set_line_b_factor($line, $b_factor);
	}
	print OUTPUT_PDB "$line\n";
}
close(OUTPUT_PDB);

# set the b-factor on each ATOM line
sub set_line_b_factor {
	my($line, $b_factor) = @_;
	my $fisrt_part = substr($line, 0, 60);
	my $second_part = sprintf("%6.2f", $b_factor);
	my $third_part = substr($line, 66);
	return $fisrt_part.$second_part.$third_part;
}

# Usage information of this script
sub Usage {
	print "Usage: set_pdb_b_factor.pl 
		-b b-factor file
		   user-provided b_factor file
		-p input pdb file
		   the pdb file in which the b-factor column will be set
		-o output pdb file
		   file in which the reset pdb will be stored
		-h help information [optional]
		   prints help information\n";
}

# Parse command line input
sub ParseInput {
	# prints out usage information if "-h" is given regardless of 
	# the fact whether others flags are specified
	for my $i (0..$#ARGV) {
		if($ARGV[$i] eq "-h") {
			Usage();
			exit;
		}
	}
	# parse in the command line arguments
	for my $j (0..$#ARGV) {
		if($ARGV[$j] eq "-b") {
			$b_factor_file = $ARGV[$j+1];
		} elsif($ARGV[$j] eq "-p") {
			$input_pdb=$ARGV[$j+1];
		} elsif($ARGV[$j] eq "-o") {
			$output_pdb=$ARGV[$j+1];
		}
	}
}
