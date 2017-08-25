#!/usr/bin/perl

use strict;
use warnings;

# load required libraries
use Getopt::Long;
use File::Basename;

# rename command line variable
my $query;
my $dbpath ||= "/dors/meilerlab/apps/PISCES/opm/uniprot20_2016_02/uniprot20_2016_02";

# get command line arguments
GetOptions(
	"query|q=s" => \$query,
	"db|d=s"  => \$dbpath,
) or die "Could not process command line options: $!";

# get PDB ID from sequence file name
my $pdb_id  = (fileparse($query, ".fasta"))[0];

# run hhblits
unless ( -s "$pdb_id.a3m" && -s "$pdb_id.psi" ) {
	my $hhblits_command =
	    "hhblits "
	  . "-i $query "
	  . "-d $dbpath " . "-n 2 "
	  . "-oa3m $pdb_id.a3m "
	  . "-o $pdb_id.hhr "
	  . "-cov 75 "
	  . "-qid 25 "
	  . "-id 90 "
	  . "-cpu 10 "
	  . "-opsi $pdb_id.psi "
	  . "-neffmax 5 "
	  . "-e 0.00001";
	print "\nRunning $hhblits_command ...\n";
	system($hhblits_command) == 0
	  or die "hhblits failed against $pdb_id...\n";
}

# reformat multiple sequence alignment
my @wc = split( /\s+/, qx/wc -l $pdb_id.psi/ );
my $num_hits = $wc[0];
if ( $num_hits > 300 ) {
	warn "Warning: too many sequences in $pdb_id.psi, "
	  . "rate4site is probably not able to handle."
	  . "Trim down to 300 sequences ...";

	# trim down to 300 sequences
	system("cp $pdb_id.psi $pdb_id.psi.old");
	system("sed -i -n  '1,300p' $pdb_id.psi");
	system("reformat.pl psi fas $pdb_id.psi $pdb_id.msa") == 0
		or die "make sure that reformat.pl can be found in your search PATH";
}
unless ( -s "$pdb_id.msa" ) {
	my $reformat2fas = "reformat.pl a3m fas $pdb_id.a3m $pdb_id.msa";
	print "Running $reformat2fas ...\n";
	system($reformat2fas) == 0
	  or die "reformat.pl failed against $pdb_id...\n";
}
unless ( -s "$pdb_id.clustal" ) {
	my $reformat2clu = "reformat.pl a3m clu $pdb_id.a3m $pdb_id.clustal";
	print "Running $reformat2clu ...\n";
	system($reformat2clu) == 0 or die "reformat.pl failed against $pdb_id...\n";
}

# run al2co
# unless ( -s "$pdb_id.al2co" ) {
#	my $al2co = "al2co -i $pdb_id.clustal -o $pdb_id.al2co -a T -g 0.5";
#	print "Running $al2co ...\n";
#	system($al2co) == 0 or die "al2co failed against $pdb_id...\n";
# }

# run rate4site
unless ( -s "$pdb_id.res" ) {
	my $rate4site =
	    "rate4site "
	  . "-s $pdb_id.msa "
	  . "-o $pdb_id.res "
	  . "-x $pdb_id.tree "
	  . "-y $pdb_id.orig.res";
	print "Running $rate4site ...\n";
	system($rate4site) == 0
	  or die "rate4site failed against $pdb_id: $!";
}

=head1 NAME

<application name> - <one line description of application's purpose>

=head1 VERSION

This documentation refers to <application name> version 0.0.1.

=head1 DESCRIPTION

A full description of the application and its features.

=head1 USAGE

	# Brief working invocation example(s) here showing the most common usage(s)
	
	# This section will be as far as many users ever read,
	# so make it as educational and exemplary as possible.

=head1 REQUIRED ARGUMENTS

A complete list of every argument that must appear on the command line.
When the application is invoked, explaining what each of them does, any
restrictions on where each one may appear, and how the various arguments 
and options may interact.

=head1 OPTIONS

A complete list of every available option with which the application
can be invoked, explaining what each does, and listing any restrictions,
or interactions.

=head1 DIAGNOSTICS

A list of very error and warning message that the application can generate, with a
full explanation of each problem, one or more likely causes, and any suggested
solutions.

=head1 AUTHOR

<Authos name(s)> (<contact address>)

=head1 COPYRIGHT

Copyright (c) <year> <copyright holder>. All rights reserved.

=head1 LICENSE

This application is free software; you can redistribute it and/or modify it under
the same terms as Perl itself.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

