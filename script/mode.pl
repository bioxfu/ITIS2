#sub mode {
#	my @nums = @_;
#	my %ha;
#	foreach my $i (@nums){
#		$ha{$i} ++;
#	}
#	my ($mode) = (sort {$ha{$a} <=> $ha{$b}} keys %ha)[-1];
#	return $mode;
#}

sub mode {
	my @nums = @_;
	my %ha;
	foreach my $i (@nums){
		$ha{$i} ++;
	}
	my $max_count = (sort{$a<=>$b} values %ha)[-1];
	print "$max_count\n";
	my @modes;
	foreach my $n (@nums) {
		print "$n\n";
		if ($ha{$n} == $max_count) {
			push @modes, $n;		
		}
	}
	my $mode = (sort @modes)[0]; 
	return $mode;
}



my @num = (1220051,1220059,1220063,1220068,1220063,1220068);

my $m = mode(@num);

print "mode = $m\n";

