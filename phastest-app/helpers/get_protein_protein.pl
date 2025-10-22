#!/usr/bin/perl -w


# open the file 'true_prophage.txt' first.
# then get the NC_XXX and region number , go into that NC_XXX folder and 
# open extract_result.txt.tmp to get the PHAGE proteins and Evalues.
use Cwd;
my $sub_program_dir = "/var/www/html/phast/current/public/sub_programs";
my $cur_dir = getcwd;
my $database_file = "$sub_program_dir/blast/data/phage_database",
my $header;

open(IN, 'true_prophage.txt') or die "Cannot open file 'true_prophage.txt'";
open(OUT, '>true_prophage_phage_protein_ORF_position_compare.txt') or die "Cannot write file 'true_prophage_phage_protein_ORF_position_compare.txt'";
while (<IN>){
	if ($_=~/^CHROMOSOME/ or $_=~/^--------------/){
		$header .= $_;
	}
	if ($_=~/^NC_/){
		print OUT $header;
		print OUT $_;
		@array=();
		@array =split(/\s\s\s+/, $_);
		$NC = $array[0];
		$region_num = $array[1];
		$most_common_phage= $array[$#array-2];
		
		chdir $cur_dir;
		#chdir  $NC;
		
		if (!(-e "$NC.txt") ){
			next;
		}
		open (IN2, "$NC.txt") ;
		%hash=();
		$start=0;
		while(<IN2>){
			if ($_ =~/region (\d+) is from/){
				$start=1;
			}
			if ($start==1){
				$hash{$1} .= $_;
			}
				
		}
		close IN2;


		open(IN3, $database_file) or die "Cannot open $database_file" ;
		%hash1=();
		my $count1 =0;
		while(<IN3>){
			if ($_=~/$most_common_phage.*gi:(\d+)/){
				$index = $1;
				#print "jjjj\n";
				$hash1{$index} = ++$count1;
				if ($_=~/(complement\(\d+\.\.\d+\))/){
					$hash1{$index} .= "\|$1";
				}elsif ($_=~/(\d+\.\.\d+)/){
					$hash1{$index} .= "\|$1";
				}
				
			}
		}
		close IN3;
		
		my $data1 ='';
		$data1 =`cat extract_result.txt`;
		my @data = split("\n", $data1);

		@lines =  split ("\n", $hash{$region_num});
		$header1 = sprintf ("\n%-40s     %-35s     %-25s     %-40s     %-25s     %-25s     %-10s\n", '', "PROPHAGE_END5_POSITION", "PROPHAGE_ORF_POSITION","BLAST_HIT_PHAGE_PROTEIN","PHAGE_ORF_POSITION", , "PHAGE_END5_POSITION","E-VALUE");
		$dash_line='';
		for($i =0; $i<=185; $i++){
			$dash_line .='-';
		}
		$dash_line = sprintf("%-40s     %s\n", '', $dash_line);
		print OUT $header1;
		print OUT $dash_line;
		my $count = 0;
		foreach $line(@lines){
			#1565207  53718978           gi|38707949|ref|NP_945090.1|, TAG = PHAGE_burkho_phi1026b, E-VALUE = 3e-68
			if ($line =~/^(\d+)/) {
				$count++;
			}
			my $a='';
			my $b='';
			my $c='';
			my $d='';
			my $e='';
			if ($line =~/^(\d+).*(PHAGE_.*),.*E-VALUE = (\S+)/){
					$a= $1; $b= $2; $c = $3;
					#print "XXXX  ".$line;
					#print "bbbbbbb=$b\n";
					@temp_array=();
					foreach $m (@data){
						if ($m=~/$a/ && $m =~/$b/){
							@temp_array = split(/\s\s\s+/, $m);
							$a = $temp_array[0];
							last;
						}
					}
					if($line =~/gi\|(\d+)\|/){
						$gi=$1;
						
						if ($b=~/$most_common_phage/ &&  length($b) == length($most_common_phage)){
							if (!defined($hash1{$gi})){
								$hash1{$gi} = 'N/A';
								
							}else{
								$d = $hash1{$gi};
								$d =~s/\|.*//;
								$e = $hash1{$gi};
								$e =~s/.*\|//;
							}
							$line = sprintf ("%-40s     %-35s     %-25s     %-40s     %-25s     %-25s     %-10s\n", '', $a, $count, "$b, gi$gi", $d, $e ,$c);
							
						}else{
							$line = sprintf ("%-40s     %-35s     %-25s     %-40s     %-25s     %-25s     %-10s\n", '', $a, '', "$b, gi$gi", '','', $c);
							
						}
					}else{
						if ($b=~/$most_common_phage/){
							$line = sprintf ("%-40s     %-35s     %-25s     %-40s     %-25s     %-25s     %-10s\n", '', $a, $count, $b, '','', $c);
						}else{
							$line = sprintf ("%-40s     %-35s     %-25s     %-40s     %-25s     %-25s     %-10s\n", '', $a, '', $b, '', '', $c);
						}
						
					}
					
					print OUT $line;
					
				
			}
		}
		print OUT "\n\n\n\n";

		
	}
}
close IN;
close OUT;
exit;

