#!/usr/bin/perl -w


# this program is used to extract true prophage regions from file 'true_defective_prophage.txt'
#intput :true_defective_prophage.txt;
# outpuot: true_prophage.txt.
open(IN, 'true_defective_prophage.txt') or die "Cannot open true_defective_prophage.txt";
open(OUT, '>true_prophage.txt') or die 'Cannot write true_prophage.txt';

my $data='';
my $start=0;

my $NC= '';
my %hash =();
my @array=();
my @temp=();
my $header ='';
my $separate_line='';
my $count=0;

while(<IN>){
	if ($_=~/\[(.*)\]/){
		
		$NC=$1;
		next;
	}
	if ($_=~/^\s+REGION/){
		$header = $_;
		$start = 1;
		next;
	}
	if ($start==1 && $_=~/^\s*$/){
		$start=0;
	}
	if ($start ==1){
		if ($_=~/---------/){
			$separate_line = $_;
		}
		if ($_=~/incomplete/){
			next;
		}
		if ($_=~/true/){
			@array =split(/\s\s\s+/, $_);
			push @temp, $array[$#array-2];
			$count++;
			$hash{$NC."\|$count"} = $_;
		}
			
	}
}
print "true_count= $count\n";
my @ar=();
my $flag =0;
foreach $a (@temp){
	$flag = 0;
	foreach $b(@ar){
		if ($a eq $b){
			$flag =1;
			last;
		}
	}
	if ($flag ==0){
		push @ar, $a;
	}
}
print OUT "CHROMOSOME".$header;
$separate_line =~s/^\s*//;
print OUT "--------------------------------".$separate_line;
my $count2=0;
foreach $a (@ar){
	foreach $k (keys %hash){
		$NC = $k;
		$NC=~s/\|.*//;
			
		if ($hash{$k} =~/$a/){
			$count2++;
			if ($a !~/transposase/ && $a!~/recombinase/ && $a!~/integrase/){ 
				print OUT  "$NC     $hash{$k}";
			}
			$hash{$k}='';
		}
	}
	print OUT "\n";
}
print "count2=$count2\n";
close IN;
close OUT;
exit;

