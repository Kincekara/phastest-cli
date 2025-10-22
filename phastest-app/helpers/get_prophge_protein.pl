#!/usr/bin/perl -w

# extract protein names with "prophage protein" from 54 NC 

my $tmp="/var/www/html/phast/current/public/tmp";

open(IN, "NC_list") or die "Cannot open NC_list";
open (OUT, ">propahge.db") or die "Cannot write prophage.db";
while(<IN>) {
	my $NC=$_;
	chomp($NC);
	chdir "$tmp/$NC";
	open (IN1, "$NC.faa") or die "Cannot open $NC.faa";
	my $flag = 0;
	while (<IN1>){
		if ($_=~/>/ ){
			if($_=~/prophage [protein|integrase|capsid|fiber|tail|plate|transposase|coat|head|portal|terminase|protease|lysis|lysin|repressor].*\[(.*)\]/i or $_=~/phage-like.*[protein|integrase|capsid|fiber|tail|plate|transposase|coat|head|portal|terminase|protease|lysis|lysin|repressor]\[(.*)\]/i) {
				my $name = $1;
				my @arr = split(" ", $name);
				my $prefix = "PROPHAGE_".substr($arr[0], 0, 6)."_$arr[$#arr]";
				$_=~s/>gi/>$prefix-gi/;
				print OUT $_;
				$flag = 1;
			}else{
				$flag = 0;
			}
			next;
		}
		if ($flag ==1){
			print OUT $_;
		}
	}
	close IN1;

}			
close IN;

exit;

