#!/usr/bin/perl -w

#get integrase from ncbi.out
my $sub_program_dir = "/var/www/html/phast/current/public/sub_programs";
my $data = `cat $sub_program_dir/phage_finder/DB/virus.db`;
open (OUT, ">combined.hmm_FRAG") or die "Cannot write combined.hmm_FRAG";
print OUT "Query HMM:   Integrase\n";
print OUT "Accession:   Integrase\n";
print OUT "Description: Integrase\n";
print OUT "Sequence                     Description                Score    E-value  N\n";
print OUT "--------                     -----------                -----    ------- ---\n";
open(IN, "ncbi.out") or die "Cannot open ncbi.out";
while(<IN>){
	my @arr = split("\t", $_);
	my $gi='';
	if ($arr[1]=~/(gi\|\d+\|)/){
		$gi=$1;
		$gi =~s/\|/\\\|/g;
	}
	my $pro_name ='';
	if ($gi ne ''){
		if ($data=~/($gi\w+\|\S+?\|.*?)\n/){
			$pro_name = $1;
		}
	}else{
		print "gi='$gi', not found in virus.db\n" ;
	}
	if ($pro_name=~/Integrase/i or $pro_name =~/Int/i ){
		$pro_name = substr($pro_name, 0, 53);
		print OUT "$pro_name   313.6    4.3e-92   1\n";
	}
}
close IN; 
print OUT "\n";
close OUT;
exit;


