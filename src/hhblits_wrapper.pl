#!/usr/bin/perl

use strict;
use warnings;

# load required libraries
use Getopt::Long qw(GetOptions);
use File::Basename qw(fileparse);

# rename command line variable
my $query;
my $dbpath ||= "/dors/meilerlab/apps/PISCES/opm/uniprot20_2016_02/uniprot20_2016_02";
my $evalue ||= 0.001;
my $num_iterations ||= 2;
chomp( my $nproc ||= `nproc` );

# get command line arguments
GetOptions(
	"query|q=s"          => \$query,
	"db|d=s"             => \$dbpath,
	"evalue|e=f"         => \$evalue,
	"num_iterations|i=i" => \$num_iterations,
	"nproc|p=i"          => \$nproc,
) or die "Could not process command line options: $!";

# get PDB ID from sequence file name
my $pdb_id = ( fileparse( $query, ".fasta" ) )[0];

# run hhblits
unless ( -s "$pdb_id.a3m" && -s "$pdb_id.psi" ) {
	my $hhblits_command =
	    "hhblits "
	  . "-i $query "
	  . "-d $dbpath "
	  . "-n $num_iterations "
	  . "-oa3m $pdb_id.a3m "
	  . "-o $pdb_id.hhr "
	  . "-cpu $nproc "
	  . "-opsi $pdb_id.psi "
	  . "-e $evalue";
	print "\nRunning $hhblits_command ...\n";
	system($hhblits_command) == 0
	  or die "hhblits failed against $pdb_id ...\n";
}

# reformat multiple sequence alignment
unless ( -s "$pdb_id.mfa" ) {
	my $reformat2fas = "reformat.pl a3m fas $pdb_id.a3m $pdb_id.mfa";
	print "Running $reformat2fas ...\n";
	system($reformat2fas) == 0
	  or die "reformat.pl failed against $pdb_id...\n";
}
unless ( -s "$pdb_id.aln" ) {
	my $reformat2clu = "reformat.pl a3m clu $pdb_id.a3m $pdb_id.aln";
	print "Running $reformat2clu ...\n";
	system($reformat2clu) == 0 or die "reformat.pl failed against $pdb_id...\n";
}
