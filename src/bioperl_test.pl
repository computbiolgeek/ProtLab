#!/bin/perl

use strict;
use warnings;

use File::Basename;
use Bio::SeqIO;



# get the filename somehow

my $file = shift;

my $seqio_obj = Bio::SeqIO->new(-file => $file);