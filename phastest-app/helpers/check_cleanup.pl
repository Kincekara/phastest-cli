#!/usr/bin/perl -w
# This program is used to check disk space and release space 
# if it is empty.
my $exec_dir = "/var/www/html/phast/current/public/cgi-bin";
my $no_space_cleanup_days = 2;
my $tmp_dir = "/var/www/html/phast/current/public/tmp";
my $space = `df -h ~/project/tmp`; 
my ($available_space) = $space =~ /(?:.*?\n)\S+\s+\S+\s+\S+\s+(\S+)/s;
my $mseg =  "Available space of ~/project/tmp is $available_space\n";
print $mseg;
if ($available_space =~/^\d+$/ && $available_space == 0) {
   print "No space on tmp folder. call cleanup.pl\n";
   my $cm1 = "perl $exec_dir/cleanup.pl  $tmp_dir $no_space_cleanup_days";
   print "$cm1\n";
   system($cm1);
}
$space = `df -h ~/project/tmp`;
($available_space) = $space =~ /(?:.*?\n)\S+\s+\S+\s+\S+\s+(\S+)/s;
print "Now space is $available_space\n";
print "Program exit\n";
exit;
