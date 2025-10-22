#!/usr/bin/perl 

# this program will transfer the result to botha-w6
#
my $NC= $ARGV[0];

chdir "/home/prion/phage/result_tmp";
if( -d $NC){
	system("scp -i ~/.ssh/scp-key -r $NC prion\@botha-w6:/local/data/phage/result_tmp/.");
}
exit;




