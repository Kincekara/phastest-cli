#!/usr/bin/perl -w

# this program is used to make genomic comparisons to each group that 
# has related phage and prophages by using clustalw.
# input: true_prophage_phage_protein_ORF_position_compare.txt
# output : all DAN files in  folder 'genomic_comparison'.

# tmp files are in folder 'genomic_comparison';
use Cwd;
use Bio::Perl;
use Bio::DB::GenBank;


my $cur_dir= getcwd;
my $tmp_dir = $cur_dir."/genomic_comparison";
my $phage_DNA_file = '/var/www/bateria_phage/PHAGE_ALL_DNA/phage_list.clw';

open(IN, 'true_prophage_phage_protein_ORF_position_compare.txt') or die "Cannot open 'true_prophage_phage_protein_ORF_position_compare.txt'";
$start =0;
my @temp=();
while(<IN>){
	if ($start ==0 && $_=~/^\s*$/){
		next;
	} 
	if ($_=~/^NC_/){
		@array = split (/\s\s\s+/, $_);
		$NC = $array[0];
		$region_num = $array[1];
		$phage = $array[$#array-2];
		$region = $array[3];
		$start =0;
		
		next;
	}
	if ($_=~/PROPHAGE_END5_POSITION                  PROPHAGE_ORF_POSITION/) {
		$start = 1;
		@temp=();
		next;
	}
	if ($start ==1 && $_=~/\.\..*\.\./ ){
		@array = $_ =~/(complement)/g;
		
		if ($#array ==-1 or $#array==1) {
			$complement = 0;
		}elsif($#array==0){
			$complement = 1;
		}
		#print "complemnt= '$complement' '$#array'  $_";
		push @temp, $complement;
		next;
	}
	if ($start ==1 && $_=~/^\s*$/){
		$start = 0;
		$count =0;
		$count1 =0;
		foreach $a (@temp){
			if ($a ==1){
				$count++;
			}else{
				$count1++;
			}
		}
		if ($count > $count1){
			$complement =1;
		}else{
			$complement = 0;
		}
		print "$NC  $region_num $region  $phage '$count' '$count1' '$complement'\n";
		#exit ;

		# handle this case;
		get_DNA_seq($NC, $region, $phage, $complement, $region_num);
		next;
	}
}
close IN;
exit;

sub get_DNA_seq {
	my ($NC, $region, $phage, $complement, $region_num) =@_;
	$phage =~s/\s/_/g;
	if (!(-d "$tmp_dir/$phage")){
		system("mkdir $tmp_dir/$phage");
	}
	if (-e "$tmp_dir/$phage/$phage\_DNA"){
		system("rm -rf $tmp_dir/$phage/$phage\_DNA");
	}
	open (IN1, $phage_DNA_file) or die "Cannot open $phage_DNA_file";
	open (OUT, ">$tmp_dir/$phage/$phage\_DNA") or die "Cannot write $tmp_dir/$phage/$phage\_DNA";
		$start =0;
		while(<IN1>){
			if ($_=~/>$phage/){
				print OUT $_;
				$start=1;
				next;
			}
			if ($start==1){
				print OUT $_;
				last;
			}
	}
	close IN1;
	close OUT;
	

	#get back the prophage ncbi file
	my $db = Bio::DB::GenBank->new();
	my $seq_object;
	if ($NC ne 'NC_001317'){
		$seq_object = $db->get_Seq_by_acc($NC);
	}else{
		$seq_object = $db->get_Seq_by_gi(9634055);
	}		
	$seq = $seq_object->seq;
	$region=~/(\d+)-(\d+)/;
	$start = $1; $end = $2;
	$seq = substr($seq, $start-1, $end-$start+1);
	my $com='';
	if ($complement==1){
		$seq = complement($seq);
		$seq = scalar reverse $seq;
		$com = 'complement';
	}
	open(OUT , ">$tmp_dir/$phage/$NC\_$region_num\_DNA") or die "Cannot write $tmp_dir/$phage/$NC\_$region_num\_DNA";
	print OUT ">$NC\_$region_num ; $region; $com\n";
	print OUT $seq;
	close OUT;


}

sub complement{
	my $str = shift;
	$str =~tr/[A, C, G, T]/[T, G, C, A]/;
	return $str;

}


