#!/usr/bin/perl 
#
#reorder the ncbi.out 
#@input: ncbi.out 
#@output: updated ncbi.out
#Usage: perl reroder_ncbi_out.pl  <ncbi.out's path>
my $file = $ARGV[0];
open(IN, $file) or die "Cannot open $file";

my %hash = ();
while(<IN>){
	#gi|00001|ref|NC_000000| PHAGE_Entero_JS98_NC_010105-gi|161622590|ref|YP_001595284.1|    91.96   286     23      0
	#       1       286     1       286     4e-149   525
	if ($_=~/^gi\|(\d+)\|/){
		$hash{$1} .= $_; 	
	}
}	
close IN;
my @array = sort(keys %hash);
open(OUT, "> $file") or die "Cannot write $file";
foreach my $k (@array){
#	print "$k\n";
	print OUT $hash{$k};
}
close OUT;
print "$file is order\nProgram exit\n";
exit(0);
