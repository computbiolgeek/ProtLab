#!/usr/bin/perl

use strict;
use warnings;

use constant GP => -2;

my ( $seq_a, $seq_b ) = qw(HEAGAWGHEE PAWHEAE);

my $alignment = get_alignment( $seq_a, $seq_b );

print ${$alignment}[0] . "\n";
print ${$alignment}[1] . "\n";
print "Score: ${$alignment}[2]\n";

######################INTERFACE SUBROUTINE##########################
# Comments     : Compute the maximum of a list of numbers.
# Usage        : max(@numbers)
# Parameter(s) : A list of numbers stored in an array.
# Returns      : The maximum of the list.
####################################################################
sub max {
	my ( $max, @numbers ) = @_; # $max gets assigned the first number in @_
	foreach my $number (@numbers) {
		$max = $number if $number > $max;
	}
	return $max;
}

######################INTERFACE SUBROUTINE##########################
# Comments     : Compute the minimum of a list of numbers.
# Usage        : min(@numbers)
# Parameter(s) : A list of numbers stored in an array.
# Returns      : The minimum of the list.
####################################################################
sub min {
	my ( $min, @numbers ) = @_; # $min gets assigned the first number in @_
	foreach my $number (@numbers) {
		$min = $number if $number < $min;
	}
	return $min;
}

######################INTERFACE SUBROUTINE##########################
# Comments     : Returns the score of a pair of residues.
# Usage        : score_aa_pair($aa_a, $aa_b)
# Parameter(s) :
#                $aa_a: first amino acid residue
#                $aa_b: second amino acid residue
# Returns      : The score the passed amino acid pair.
####################################################################
sub score_aa_pair {
	my ( $aa_a, $aa_b ) = @_;
	return ( $aa_a eq $aa_b ) ? 2 : -1;
}

######################INTERFACE SUBROUTINE##########################
# Comments     : Computes the alignment score of two sequences using
# 			     the Needleman-Wunsch algorithm.
# Usage        : score_alignment($sequence_a, $sequence_b)
# Parameter(s) :
#                $sequence_a: horizontal sequence
#                $sequence_b: vertical sequence
# Returns      : The optimal alignment score of two sequences.
####################################################################
sub score_alignment {
	my ( $sequence_a, $sequence_b ) = @_;

	# matrix for storing dynamic programming scores
	my @dp_matrix = ();

	# populate the first row and first column of the matrix
	foreach my $i ( 0 .. length($sequence_a) ) {
		$dp_matrix[0][$i] = GP * $i;
	}
	foreach my $j ( 0 .. length($sequence_b) ) {
		$dp_matrix[$j][0] = GP * $j;
	}

	# compute score for each entry according to Needleman-Wunsch
	foreach my $i ( 1 .. length($sequence_a) ) {
		foreach my $j ( 1 .. length($sequence_b) ) {
			my $pair_score = score_aa_pair(
				substr( $sequence_a, $i - 1, 1 ),
				substr( $sequence_b, $j - 1, 1 )
			);
			$dp_matrix[$j][$i] = max(
				$dp_matrix[ $j - 1 ][$i] + GP,
				$dp_matrix[$j][ $i - 1 ] + GP,
				$dp_matrix[ $j - 1 ][ $i - 1 ] + $pair_score
			);
		}
	}

	return \@dp_matrix;
}

######################INTERFACE SUBROUTINE##########################
# Comments     :
# Usage        :
# Parameter(s) :
# Returns      :
####################################################################
sub get_alignment {
	my ( $sequence_a, $sequence_b ) = @_;

	# compute the scores
	my $dp_matrix = score_alignment( $sequence_a, $sequence_b );

	my ( $i, $j ) = ( length($sequence_a), length($sequence_b) );

	# best score
	my $best_score = ${$dp_matrix}[$j][$i];

	# array that stores the two aligned sequences
	my ( $aligned_a, $aligned_b ) = ( "", "" );
	while ( $i > 0 && $j > 0 ) {

		# if the current score comes from the diagonal, align the two residues
		if (
			${$dp_matrix}[$j][$i] ==
			${$dp_matrix}[ $j - 1 ][ $i - 1 ] + score_aa_pair(
				substr( $sequence_a, $i - 1, 1 ),
				substr( $sequence_b, $j - 1, 1 )
			)
		  )
		{
			$aligned_a .= substr( $sequence_a, $i - 1, 1 );
			$aligned_b .= substr( $sequence_b, $j - 1, 1 );
			--$i;
			--$j;
		}

		# if the score comes vertically, add a gap to $aligned_a
		elsif ( ${$dp_matrix}[$j][$i] == ${$dp_matrix}[ $j - 1 ][$i] + GP ) {
			$aligned_a .= "-";
			$aligned_b .= substr( $sequence_b, $j - 1, 1 );
			--$j;
		}

		# if the score comes horizontally, add a gap to $aligned_b
		else {
			$aligned_a .= substr( $sequence_a, $i - 1, 1 );
			--$i;
			$aligned_b .= "-";
		}
	}

	# if either $sequence_a or $sequence_b has residues left, align them to gaps
	while ( $i > 0 ) {
		$aligned_a .= substr( $sequence_a, $i - 1, 1 );
		--$i;
		$aligned_b .= "-";
	}
	while ( $j > 0 ) {
		$aligned_a .= "-";
		$aligned_b .= substr( $sequence_b, $j - 1, 1 );
		--$j;
	}
	$aligned_a = reverse($aligned_a);
	$aligned_b = reverse($aligned_b);
	return [ $aligned_a, $aligned_b, $best_score ];
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