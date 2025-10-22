#!/usr/bin/perl -w



open(IN, $ARGV[0]) or die "Cannot open $ARGV[0]";
open(OUT, ">$ARGV[0].tmp" ) or die "Cannot write $ARGV[0].tmp";
my %hash =();
my %hash1=();
my $flag = 0;
while(<IN>) {
	if ($_=~/>gi\|(\d+)/){
		#>gi|158333234|ref|YP_001514406.1| NUDIX hydrolase [Aca
		$hash{$1} += 1;
		$hash1{$1} .= "$_\n";
		if ($hash{$1}>1){
			$flag = 1;
		}else{
			$flag =0;
			print OUT  $_;
		}

	}else{
		if ($flag==0){
			print OUT $_;
		}
	}
		
}
close IN;
close OUT;
my $count = 0;
foreach my $k (keys %hash){
	if ($hash{$k}>1){
		print $hash1{$k};
		$count++;
	}
}
print "count=$count\n";
exit;

