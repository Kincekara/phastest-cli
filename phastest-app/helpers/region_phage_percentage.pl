#!/usr/bin/perl -w

# this program will get the percentage of each phage proteins against the protein of that prophage region
# extract_result.txt needed!
#input: list file of prophage folder
#output: region_phage_percentage.txt
# Usage : perl region_phage_percentage.pl <NC_number> 
use Cwd;
use Bio::Perl;
use Bio::DB::GenBank;
use Data::Dumper;

open(OUT, '>region_phage_percentage.txt') or die "Cannot write region_phage_percentage.txt";
get_percentage(\*OUT, "$ARGV[0]\_phmedio.txt");
close OUT;
open(IN, "region_phage_percentage.txt") or die "Cannot open region_phage_percentage.txt";

my %hhash=();
my $header='';
my $head_flag = 1;
my $p='';
my $region='';
while(<IN>){
	if ($_=~/\s+REGION\s+POSITION/){
		if ($p ne ''){
			$hhash{$region}= $p;
		}
		$p=$_;
		$head_flag=0;
	}
	if ($head_flag==1){
		$header.=$_;
		next;
	}
	if ($_=~/-------------------/){
		$p.=$_;
		next;
	}
	if ($_=~/\s+(\d+)\s+/){
		$region=$1;
		$p.=$_;
	}
}
close IN;
if ($p ne ''){
	$hhash{$region}= $p;
} 
open(OUT, '>region_phage_percentage.txt.tmp');
print OUT $header;
foreach my $key (sort {$a<=>$b} keys %hhash){
	print "XXX $key\n";
	print OUT $hhash{$key};
}	
close OUT;
system("mv region_phage_percentage.txt.tmp  region_phage_percentage.txt");
exit;

sub get_percentage{
	my $fh = shift;
	my $scan_output_file = shift;
	my $flag = -1;
	my %hash1=();
	my %prot_counts=();
	if (-e $scan_output_file){
		get_perc($scan_output_file, $fh, \$flag, \%hash1, \%prot_counts);
		#print Dumper(\%hash1);
	}
	$region ='';
	$end5_3 = '';
	$prot_name='';
	$perc = '';
	$last_region='';
	
	my %hash_partial=();
		
		foreach $k (sort {$a cmp $b} keys %hash1){
			
				$region= $k ;
				$region=~s/\|(.*)//;
				$end5_3 = $1;
				
				if ($last_region ne $region){
					if (scalar(keys %hash_partial) >0){
						#print Dumper(\%hash_partial);
						foreach my $key (sort {$a<=>$b} keys  %hash_partial){
							print $fh $hash_partial{$key};
						}
					}
					$line = sprintf ("\n%10s %-8s   %-20s   %-100s   %-10s   %-10s   %-10s\n", "", "REGION", "POSITION", "PHAGE", "NUM", "TOTAL_GENE_NUM", "PERCENTAGE");
					print $fh $line;
					$line =sprintf("%10s %s\n", "",  "-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------");
					print $fh $line;
					$last_region = $region;
					%hash_partial=();
				}	
			 	$begin = $end5_3;
				$begin =~s/_\d+//;
				$prot_name = $hash1{$k};
				$prot_name =~s/\|(\d+)$//;
				my $test=$1;
				if ($test !~/^[\d.]+$/ or $prot_counts{$region} !~/^[\d.]+$/){
					#print "DDDD hash{k}='$hash1{$k}', test = '$test' , prot_counts{region} =  $prot_counts{$region}\n";
					#exit;
				}
				$perc = int($test/$prot_counts{$region}*1000);
				$perc = $perc/10;
				$num = $test;
				if ($prot_name!~/attachment_site/){
					$line = sprintf ("%10s %-8s   %-20s   %-100s   %-10s   %-10s   %-10s\n", "", $region, $end5_3, $prot_name, $num, $prot_counts{$region}, "$perc%");
				}else{
					$line = sprintf ("%10s %-8s   %-20s   %-100s   %-10s   %-10s   %-10s\n", "", $region, $end5_3, $prot_name, '-', '0', '0');
				}
				$hash_partial{$begin} = $line;
			
		}	
	if (scalar(keys %hash_partial) >0){
		#print Dumper(\%hash_partial);
		foreach my $key (sort {$a<=>$b} keys  %hash_partial){
			print $fh $hash_partial{$key};
		}
	}
	print $fh "\n";
}

sub get_perc{
	my $f = shift;
	my $fh = shift;
	my $flag = shift;
	my $hash1=shift;
	my $prot_counts=shift;
	my %hash =();
	my $region ='';
	my $start='';
	my $end='';
	my $end5='';
	my @array =();
	my @regions=();
	my $prophage_name = '';
	my $have_HMM_flag=0;
	my $protein_name='';
	my $dat = `cat extract_result.txt`;
	open (IN, $f) or die "Cannot open $f";
	while (<IN>){
		#gi|191639869|ref|NC_011002.1| Burkholderia cenocepacia J2315 chromosome
		if ($_=~ /gi\|\d+\|ref\|\S+\|(.*),.*gc\%/ or $_=~ /gi\|\d+\|(.*),.*gc\%/ ){
			$prophage_name = $1;
			$prophage_name=~s/DEFINITION\s*//;  
			next;
		}
		if ($_=~/region (\d+) is from (\d+) to (\d+)/){
			$region =$1;
			$start = $2;
			$end = $3;
			$end5='';
			push @regions,$region;
			next;
		} 
		
		if ($_=~/^(\d+)\s+\d+/ or $_=~/^(\d+)\s+attL/ or $_=~/^(\d+)\s+attR/ or $_=~/^(\d+)\s+tRNA/ or $_=~/^(\d+)\s+tmRNA/){
			$end5 = $1;
			if ($dat =~/complement\((\d+)\.\.$end5\)/s && $end5 !=$end ){
				$end5 = $1;
			}
			elsif ($dat =~/$end5\.\.(\d+)/s  ){
				if ($1==$end){
                                $end5 = $1;
				}
                        }
			$have_HMM_flag=0;
			my $key = "$region|$end5";
			if ( $_=~/(\S*PHAGE_.*?),/  ) { #or $_=~/(\w+ phage.*);/ or $_=~/(\w* transposase.*);/ or $_=~/(\w* integrase.*);/ or $_=~/(\w* recombinase.*);/ or $_=~/\[(TRNA)\]/ or $_=~/(capsid)/i or $_=~/(tail)/i or $_=~/(fiber)/i or $_=~/(coat)/i or $_=~/(plate)/i or $_=~/(head)/i or $_=~/(prophage.*?)\[/i
				$c =$1;
				$c=~s/^\s*//;
				$protein_name ='';
				$have_HMM_flag=1;
				if ($_=~/gi\|(\d+)\|/){
					#$protein_name = get_pro_name($1);
					$protein_name = <IN>;
					if ($protein_name =~/^\s+\[ANNO\](.*)/){
						$protein_name = $1;
					}else{
						$protein_name = '';
					} 
					$hash{$key} = "$c:$protein_name";
				}else{
					$hash{$key} = $c;
				}
			}elsif ($_=~/hypothetical/){
				$hash{$key} = '     Hypothetical_protein';
			}elsif ($_=~/attR/ or $_=~/attL/){
				$hash{$key} = '     attachment_site';
			}elsif ($_=~/tRNA/ or $_=~/tmRNA/){
				$hash{$key} = '     tRNA';
			}elsif ($_=~/tmRNA/){
				$hash{$key} = '     tmRNA';
			}else{
				$hash{$key} = '     Bacterial_protein';
			}
		}
		if ($_ =~/^\s+\[HMM.*?\].*?:(.*?),/ && $have_HMM_flag==1){
			#print "$_";
			$protein_name .= $1;
			#print "   protein_name=$protein_name\n";
			$have_HMM_flag=0;
			$hash{$region."|$end5"} = $c.":$protein_name";
		}
	}
	close IN;
	my @tmp=();
	my @tmp1=();
	my @tmp2=();
	foreach(keys %hash){
		#$_ =~s/.*\|//;
		push @tmp, $_;
		#print "VVV '$_', '$hash{$_}'\n";
		
	}
	
	my $last_reg=-1;
	my @arr=();
	foreach my $key (sort {$a cmp $b}  @tmp) {
		
		@arr=split(/\|/, $key);
		if ($arr[0] != $last_reg ){
			@tmp2 = sort {$a<=>$b} @tmp2;
			foreach (@tmp2){
				$_="$last_reg|$_";
				push @tmp1, $_;
			}
			@tmp2=();
			$last_reg= $arr[0];
		}
		push @tmp2, $arr[1];
		 
	}
	
	
	if ((scalar @tmp2)>0){
		@tmp2 = sort {$a<=>$b} @tmp2;
		foreach (@tmp2){
			$_="$last_reg|$_";
			push @tmp1, $_;
		}
		@tmp2=();
		$last_reg= $arr[0];
	}
	@tmp=@tmp1;
	
	
	my $protein_count=0;
	my $phage_count=0;
	my $last_pro='N/A';
	my $last_pro_count=0;
	$start='';
	$end='';
	
	$region ='';
	$end5='';
	
	my $last_region=-1;
	my $cur_pro='';
	$protein_count=0;
	my $index = -1;
	foreach $a (@tmp) {
		$index++;
		foreach $k ( keys %hash){
			if ($k eq $a ){
				$region= $k ;
				$region=~s/\|(.*)//;
				$end5 = $1;
				if ($last_region != $region){
					if ($last_region != -1){
						$hash1->{"$last_region\|$start\_$end"}= "$last_pro\|$last_pro_count";
						#print "$last_region\|$start\_$end"."||". "$last_pro\|$last_pro_count\n";
					}
					$last_region = $region;
					$last_pro = 'N/A';
				}
				
				$prot_counts->{$region}++  if ($hash{$k}!~/attachment_site/ && $hash{$k}!~/tRNA/ && $hash{$k}!~/tmRNA/);
				if ($hash{$k} ne $last_pro && $last_pro eq 'N/A'){
					$start = $end5 ;
					$end = $end5;
					$last_pro_count=1;
					$last_pro = $hash{$k};
					#print "start=$start; k=$k; 1= '$hash{$k}'; 2= '$last_pro'\n" if ($region==1);
				}elsif ($hash{$k} ne $last_pro){
					#print "0= '$k'; 1= '$hash{$k}' ; 2= '$last_pro'\n" if ($region==1);
					$end = $end5;
					$hash1->{"$region\|$start\_$end"}= "$last_pro\|$last_pro_count";
					#print "  end=$end;  $region\|$start\_$end"."||". "$last_pro\|$last_pro_count\n" if ($region==1);
 					$last_pro = $hash{$k};
					$start = $end5 ;	
					$last_pro_count=1;
				}else{
					$end = $end5;
					$last_pro_count++;
				}
			}
		}
		if ($index == $#tmp){
		$end = $end5;
		$hash1->{"$region\|$start\_$end"}= "$last_pro\|$last_pro_count";
		#print "$region\|$start\_$end"."||". "$last_pro\|$last_pro_count\n";
		}
	}
	
=pod
	foreach $a (sort {$a cmp $b} keys %hash1) {
		print   "$a     $hash1{$a}\n";
	}
=cut
		
	print $fh  $prophage_name."\n" if ($$flag ==-1);
	$$flag = 1;
	
}

sub get_pro_name{
	my $gi_num = shift;
	#get back the prophage ncbi file
	my $db = Bio::DB::GenBank->new();
	my $seq_object = $db->get_Seq_by_gi($gi_num);
	return $seq_object->desc();
}




