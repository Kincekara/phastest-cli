#!/usr/bin/perl -w


# this program will run 54 tests cases 
# then calculate sensitivity and PPV

if (scalar @ARGV !=1){
	print STDERR "Usage: perl run_test.pl  <scan.pl>\n";
	exit(-1);
}
my $sca_exec= $ARGV[0];
print "$sca_exec\n";
chdir "~/project/tmp_test";
open (IN, "qq_test") or die "cannot open z_result.txt";
my $c = 0;
while (<IN>){
	$c++;
	chomp ($_);
	print "\t$c  $_  ".`date`."\n";
	system("perl /apps/phast/project/cgi-bin/phage_test.pl  $_  $sca_exec 2>&1|cat > /dev/null");
}
close IN;
system("perl /apps/phast/project/cgi-bin/sensitivity_ppv.pl ");
exit;
