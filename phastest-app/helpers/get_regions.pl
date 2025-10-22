#!/usr/bin/perl -w


# this program will add the regions from phpico phmedio phregions to extract_result.txt
# and  generate extract_result.txt.tmp file

use Cwd;

my $dir = $ARGV[0];

my %hash =();
if (-e "$dir\_phpico.txt"){
	open(IN2,"$dir\_phpico.txt");
	get_regions(\*IN2, \%hash);
	close IN2;
}
if (-e "$dir\_phmedio.txt"){
	open(IN3, "$dir\_phmedio.txt");
	get_regions(\*IN3, \%hash);
	close IN3;
}
if (-e "$dir\_phregions.txt"){
	open(IN4, "$dir\_phregions.txt");
	get_regions(\*IN4,\%hash);
	close IN4;
}
foreach my $k(keys %hash){
	#print "    k=$k    $hash{$k}\n";
	if ($hash{$k} =~/^(\d+).*?(\d+),$/){
		$hash{$k} = "$1,$2";
	}
	#print "    	k=$k    $hash{$k}\n";
}
open (IN5, 'extract_result.txt') or die "Cannot open extract_result.txt";
open (OUT, '>extract_result.txt.tmp') or die "Cannot write extract_result.txt.tmp";
$flag =0;
$last_region='';
$region ='';
#$first_time=1;
while(<IN5>){
	if ($_ =~/CDS_POSITION/ or $_ =~/------------------------------------------------/){
		$flag = 1;
		print OUT $_;
		next;
	}
	if ($flag ==1 ){
		my $start = '';
		if ($_=~/complement\((\d+)\.\.(\d+)\)/){
			$start=$2;
		}elsif($_=~/(\d+)\.\.(\d+)/){
			$start = $1;
		}
		
		$region =inside_region($start, \%hash);
		if ($region eq ''){
			next;
		}
		#print "   region = $region , $start=$start, last_region=$last_region\n";
		if ( $last_region ne $region ){
			print  OUT "\n#### region $region ####\n";
			$last_region = $region;
		}
	}
	if ($_=~/POSITION\s+DNA_sequence/){
		$flag =0;
	}
	print  OUT $_;
}
close IN5;
close OUT;
if (-e 'extract_result.txt.tmp'){
	change_true_defect_file();
}

exit;

# change the region position number for true_defective_prophage.txt
sub change_true_defect_file{
	my %hash = ();
	my $key ='';
	open (IN, 'extract_result.txt.tmp') or die "Cannot open extract_result.txt.tmp";
	while (<IN>) {
		if($_=~/#### region (\d+) ####/){
			$key = $1;
			next;
		}
		if($_=~/^(\d+\.\.\d+)/ or $_=~/^complement\((\d+\.\.\d+)\)/){
			$hash{$key} .= "$1,";
		}
	}
	close IN;
	foreach my $k (keys %hash){
		$hash{$k} =~s/^(\d+).*\.\.(\d+),$/$1-$2/;
	}
	open(IN, "true_defective_prophage.txt") or die "Cannot open true_defective_prophage.txt";
	open(OUT,">true_defective_prophage.txt.tmp") or die "Cannot write true_defective_prophage.txt.tmp";
	while(<IN>){
		if ($_=~/^\s+(\d+)\s+[\d\.]+Kb/){
			my $k = $1;
			$_=~ s/\d+-\d+/$hash{$k}/;
		}
		print OUT $_;
	}
	close OUT;
	close IN;
	system("mv -f true_defective_prophage.txt.tmp  true_defective_prophage.txt");
}

sub inside_region{
	my $start=shift;
	my $hash=shift;
	foreach my $k ( keys %$hash){
		my $st ='';
		my $en = '';
		($st,$en) =split("," , $hash->{$k});
		#print "      st=$st, en=$en, $start\n";
		if ($start >= $st && $start <= $en){
			if ($start == $en){
				$hash->{$k} = "-1,-1";
			}
			return $k;
		}
	}
	return '';
}


sub get_regions{
	my $fh = shift;
	my $hash=shift;
	my $region = '';
	while(<$fh>){
		#region 1 is from
		if ($_=~/region (\d+) is from/){
			$region = $1;
			next;
		}
		if ($_=~/^(\d+)/){
			#print "   region=$region   $1\n";
			$hash->{$region} .="$1,";
		}
	}
	

}


