#!/usr/bin/perl

use strict;
use warnings;

use Data::Dumper;
use Getopt::Long;

my $dbpath ||= "/dors/meilerlab/apps/PISCES/opm/2014";
my $pdbids_file;
my $log_suffix ||= ".uniref50.ascii5.log";
my $blastdb ||=
  "/dors/meilerlab/apps/BLAST/databases/db.21May2014/uniref50nofrag";
my $BLASTDBCMD = "blastdbcmd";

GetOptions(
	"path|p=s"    => \$dbpath,
	"ids|i=s"     => \$pdbids_file,
	"suffix|s=s"  => \$log_suffix,
	"blastdb|d=s" => \$blastdb,
) or die "Could not process command line options: $!";

unless ( defined $pdbids_file ) {
	die "a list of PDB IDs in a file must be supplied: $!";
}

# read in all the PDB IDs for which PSI-BLAST hits are to be extracted
open( my $fh, "<", $pdbids_file )
  or die "Could not open $pdbids_file for reading: $!";
chomp( my @pdb_ids = <$fh> );
close($fh);

# extract hits for each given PDB ID
foreach my $pdb_id (@pdb_ids) {
	my $mid_id = substr( $pdb_id, 1, 2 );
	my $psiblast_log = $dbpath . "/" . $mid_id . "/" . $pdb_id . $log_suffix;
	print "Extracting hits from $psiblast_log\n";

	# open psiblast log file to extract hits from the last iteration
	open( my $ifh, "<", $psiblast_log )
	  or die "Could not open $psiblast_log for reading: $!";
	chomp( my @lines = <$ifh> );

	# extract hits from each round
	my $hits;
	print "Start extraction from the last round, and go backwards\n";
	for ( my $i = 5 ; $i > 0 ; --$i ) {
		$hits = hits_from_round( \@lines, $i );

		# it could be that no round has more than 4 hits
		if ( scalar @{$hits} >= 4 ) {
			print "PSI-BLAST found more than 4 hits in round $i, skip other rounds\n";
			last;
		}
	}

	# warn the user if PSI-BLAST found very few hits
	if ( scalar @{$hits} < 4 ) {
		warn "Warning: PSI-BLAST found less than 4 hits in the first round, skip $pdb_id";
		next;
	}

	# open file for writing
	my $output_file = $pdb_id . "_psiblast_hits.txt";
	open( my $ofh, ">", $output_file )
	  or die "Could not open $output_file for writing: $!";
	foreach my $hit ( @{$hits} ) {
		print $ofh $hit . "\n";
	}
	close($ofh);

	# run blastdbcmd to get sequences
	my $seq_file = $pdb_id . ".seqs";
	my @command  = (
		$BLASTDBCMD,  "-db",  $blastdb,  "-entry_batch",
		$output_file, "-out", $seq_file, "-outfmt",
		"%f"
	);
	print "Running @command\n";
	system(@command) == 0 or die "blastdbcmd failed against $pdb_id: $!";
}

sub hits_from_round {
	my ( $lines, $round ) = @_;
	my @hits;
	my $is_in_round = 0;
	foreach my $line ( @{$lines} ) {
		if ( $line eq "Results from round $round" ) {
			$is_in_round = 1;
			next;
		}
		if ( !$is_in_round ) {
			next;
		}
		else {
			if ( $line =~ /(\bUniRef50_[A-Z0-9-]+\b)/ ) {
				push( @hits, $1 );
			}
			elsif ( $line =~ /^Results from round/ ) {
				last;
			}
		}
	}
	return \@hits;
}

exit;
