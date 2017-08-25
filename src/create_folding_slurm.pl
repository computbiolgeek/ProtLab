#!/usr/bin/perl

use strict;
use warnings;

use Getopt::Long qw(GetOptions);

my $pdbid;
my $dir;
my $exposure;
my $output;

GetOptions(
	"input|i=s"  => \$pdbid,
	"dir|d=s"    => \$dir,
	"exposure"   => \$exposure,
	"output|o=s" => \$output,
) or die "Could not parse command line options: $!";

unless ( defined $pdbid && defined $output ) {
	print "Usage: $0 -i <five-letter pdb id> -d <path> -o <output filename>\n";
	die "Five-letter PDB IDs for benchmark proteins must be given!";
}

my $shbang = "#!/bin/bash" . "\n";
my $job_specs =
    "#SBATCH --nodes=1" . "\n"
  . "#SBATCH --ntasks=1" . "\n"
  . "#SBATCH --mem-per-cpu=1G" . "\n"
  . "#SBATCH --time=5:00:00" . "\n"
  . "#SBATCH --array=0-99" . "\n";
my $log = "#SBATCH -o $dir/logs/%A_%a.log" . "\n";
my $lib =
    'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:'
  . '/dors/meilerlab/apps/Linux2/x86_64/lib/:'
  . '/dors/meilerlab/apps/Linux2/x86_64/lib64/' . "\n";

my $seed = 'seed = $SLURM_JOBID' . "\n";

# bcl command line options and arguments
my $bcl            = "/dors/meilerlab/apps/Linux2/x86_64/bin/bcl.exe";
my $app            = "protein:Fold";
my $fasta          = "-fasta $dir/sspred/$pdbid.fasta";
my $pool           = "-pool $dir/sspred/$pdbid.SSPredThreshold_MASP.pool";
my $native         = "-native $dir/sspred/$pdbid.pdb";
my $quality        = "-quality RMSD GDT_TS -superimpose RMSD -sspred PSIPRED MASP";
my $stage          = "-stages_read $dir/stages/stages.txt";
my $storage        = "-protein_storage $dir/$pdbid Overwrite";
my $tm_helices     = "-tm_helices $dir/sspred/$pdbid.SSPredThreshold_MASP.pool";
my $pool_prefix    = "-pool_prefix $dir/sspred/$pdbid";
my $sequence_data  = "-sequence_data $dir/sspred/ $pdbid";
my $other_stuff    = '-prefix $seed -random_seed $seed -opencl Disable';
my $score_exposure = "-score_exposure";

# bcl command line
my $command_line =
    $bcl . " " 
  . $app . " " 
  . $fasta . " " 
  . $pool . " "
  . $tm_helices . " "
  . $native . " "
  . $stage . " "
  . $quality . " "
  . $storage . " "
  . $pool_prefix . " "
  . $sequence_data . " "
  . $other_stuff;

if ($exposure) {
	$command_line .= " " . $score_exposure;
}

# output
open( my $ofh, ">", $output ) or die "Could not open $output for writing: $!";
print $ofh $shbang . $job_specs . $log . $seed . $lib . $command_line;
close($ofh);

