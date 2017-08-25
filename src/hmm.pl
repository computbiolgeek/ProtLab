#!/usr/bin/perl

use strict;
use warnings;

my @states = ( 'Rainy', 'Sunny' );

my @observations = ( 'walk', 'shop', 'clean' );

my %start_probability = (
	Rainy => 0.6,
	Sunny => 0.4,
);

my %transition_probability = (
	Rainy => { Rainy => 0.7, Sunny => 0.3, },
	Sunny => { Rainy => 0.4, Sunny => 0.6, },
);

my %emission_probability = (
	Rainy => { walk => 0.1, shop => 0.4, clean => 0.5, },
	Sunny => { walk => 0.6, shop => 0.3, clean => 0.1, },
);

######################INTERFACE SUBROUTINE##############################
# Comments     :
# Usage        :
# Parameter(s) :
# Returns      : The most likely hidden state sequence X = {x1,x2,...}
########################################################################
sub viterbi {
	
}
