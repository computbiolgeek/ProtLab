package src::math;

######################INTERFACE SUBROUTINE##########################
# Comments     : 
# Usage        : 
# Parameter(s) : 
# Returns      : 
####################################################################
sub max {
	my $max = shift @_;
	foreach (@_) {
		if ($_ > $max) {
			$max = $_
		}
	}
	return $max
}

######################INTERFACE SUBROUTINE##########################
# Comments     : 
# Usage        : 
# Parameter(s) : 
# Returns      : 
####################################################################
sub min {
	my $min = shift @_;
	foreach (@_) {
		if ($_ < $min) {
			$min = $_;
		}
	}
	return $min;
}

######################INTERFACE SUBROUTINE##########################
# Comments     : 
# Usage        : 
# Parameter(s) : 
# Returns      : 
####################################################################
sub cumsum {
	my $sum = shift @_;
	foreach (@_) {
		$sum += $_;
	}
	return $sum;
}

######################INTERFACE SUBROUTINE##########################
# Comments     : 
# Usage        : 
# Parameter(s) : 
# Returns      : 
####################################################################
sub cumprod {
	my $prod = shift @_;
	foreach (@_) {
		$prod *= $_;
	}
	return $prod;
}
