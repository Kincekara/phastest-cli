#!/usr/bin/perl -w

## this program is used to check the sensitivity and PPV
chdir "~/project/tmp_test"; 
open (IN, "qq_test_bk") or die "Cannot open qq_test_bk";
my %new_hash=();
my $inte=0;
my $att=0;
my $tRNA=0;
while(<IN>){
  my($NC)= $_=~/-a\s+(.*)/;
  if (-d $NC){
     my $data=`cat $NC/summary.txt`;
      #print "bef: ". length($data)."\n";
      $data =~s/^.*?\-\-\-\-\-\-\-\-\-\-\-\-+//s;
      #print "aft: ".length($data)."\n";
     my @arr= $data=~/\((\d+)\).*?(\d+\s*-\s*\d+)/sg;
	my @arr1=();
	 for(my $i=0; $i<scalar @arr; $i=$i+2){
		if ($arr[$i] >=0){
			push @arr1, $arr[$i+1];
		}
	 }
     $new_hash{$NC}= \@arr1;
	 my $data1=`cat $NC/detail.txt`;
	 my @att=$data1=~/(attL|attR)/gs;
	 $att += @att;
	 my @trna=$data1=~/(tRNA)/gs;
	$tRNA +=@trna;
	my @inte=$data1=~/(integrase)/gis;
	$inte += @inte;

  }
}
close IN;

open(IN, "z_result.txt") or die "Cannot open result.txt";
my $start = 0;
my %phast=();
my %standard=();
my %prophinder=();
my %phage_finder=();
my $NC='';
while (<IN>){
  next if ($_=~/^\s*$/);
  if ($_=~/(NC_\d+)/){
    $NC = $1 ;
    $phast{$NC}=[]; $standard{$NC}=[];
    $prophinder{$NC}=[]; $phage_finder{$NC}=[];
    next;
  }
  if ($_=~/reference\s+PHAST\s+Prophinder\s+Phage_finder/){
    $start =1;
    next;
  }
  next if ($start==1 && $_=~/\.\.\.\.+/);
  
  if ($start==1 && $_=~/====================+/){
    $start=0;
    next;
  }
  next if ($start==0);
  if ($start==1){
    my $s= substr($_, 0, 18); $s=~s/^\s*//; $s=~s/\s*$//;
    my $phast= substr($_, 18, 17);$phast=~s/^\s*//; $phast=~s/\s*$//;
    my $pro = substr($_, 35, 17);$pro=~s/^\s*//; $pro=~s/\s*$//;
    my $p_f = substr($_, 52);$p_f=~s/^\s*//; $p_f=~s/\s*$//;
    #print "'$s', '$phast', '$pro', '$p_f'\n";
    push @{$standard{$NC}}, $s;
    push @{$phast{$NC}}, $phast;
    push @{$prophinder{$NC}}, $pro;
    push @{$phage_finder{$NC}}, $p_f;
  }

}
close IN;

foreach my $k (keys %phast){
#    print "$k\n";
#    print "@{$phast{$k}}\n";
}

print "inte=$inte, att=$att, tRNA=$tRNA\n";

print "New phast :\n";
get_sen_ppv(\%standard, \%new_hash);
print "original phast :\n";
get_sen_ppv(\%standard, \%phast);
print "prophinder:\n";
get_sen_ppv(\%standard, \%prophinder);
print "phage_finder:\n";
get_sen_ppv(\%standard, \%phage_finder);
=pod
my $co = 0;
foreach my $k (keys %new_hash){
	 my @tmp=();
	 foreach my  $hit (@{$new_hash{$k}}){
       next if ($hit !~/\d+\s*-\s*\d+/);
       my ($start, $end)= $hit =~/(\d+)\s*-\s*(\d+)/;
		my $hit_flag = 0;
		foreach my $e (@{$phast{$k}}){		
			 next if ($e!~/\d+\s*-\s*\d+/);
	        my ($start1, $end1)= $e =~/(\d+)\s*-\s*(\d+)/;
			if ($start < $end1 && $end > $start1){
				 $hit_flag =1;
			}
		}
		if ($hit_flag ==0){
			push @tmp, $hit;
			$co++;
		}

	}
	if (scalar @tmp !=0){
		print "$k\n\t@tmp\n";
	}
}
print "count = $co\n";
=cut
exit;

sub get_sen_ppv {
  my ($standard, $query)= @_;
   my $hit_count_on_ref=0;
	my $hit_count_on_query=0;
  my $bp_hit = 0;
   my $query_bp_hit=0;
  my $query_count=0;
  my %hash=();
	my %hash1=();
  foreach my $NC (keys %$query){
     foreach my  $hit (@{$query->{$NC}}){
       next if ($hit !~/\d+\s*-\s*\d+/);
       my ($start, $end)= $hit =~/(\d+)\s*-\s*(\d+)/;
       if (!defined $start){
          print STDERR "no start defined on line '$hit, $NC'\n";
          exit;
       }
       if (! defined $end){
          print STDERR "no end defined on line '$hit, $NC'\n";
          exit;
      }
      $query_count++;
      #print "$start, $end, ";
      $query_bp_hit += $end-$start+1 if ($end-$start+1 >0);
      #print "$query_bp_hit\n";
      foreach my $e (@{$standard->{$NC}}){
        next if ($e!~/\d+\s*-\s*\d+/);
        my ($start1, $end1)= $e =~/(\d+)\s*-\s*(\d+)/;
         if (!defined $start1){
            print STDERR "no start defined in reference on  line '$e, $NC'\n";
            exit;
         }
         if (! defined $end1){
            print STDERR "no end defined in reference on line '$e, $NC'\n";
            exit;
          }
        #print "$start , $end \\ $start1, $end1, $NC\n";
        if ($start < $end1 && $end > $start1 ){
			if (!defined $hash1{$e}){
				$hash1{$e}=1;
			}else{
				$hash1{$e}++;
			}
			$hit_count_on_query++;
          my $start2= ($start1 > $start)? $start1 : $start;
          my $end2= ($end1 > $end)?  $end : $end1;
          $bp_hit += $end2 -$start2+1;
        #print "hitt\n";
          $hash{$start1}="$start , $end \\ $start1, $end1, $NC\n";
          last;
        }
      }
    }
  }
	$hit_count_on_ref = scalar keys %hash1;
  foreach my $k (sort keys %hash){
#    print $hash{$k};
  }
  my $ref_count=0;
  my $ref_bp_count = 0;
  foreach my $NC (%$standard){
      foreach my $e (@{$standard->{$NC}}){
        next if ($e!~/\d+\s*-\s*\d+/);
        my ($start, $end)= $e =~/(\d+)\s*-\s*(\d+)/;
        $ref_count++;
        $ref_bp_count +=  $end-$start+1;
      }
  }
 
  my $sen = $hit_count_on_ref/$ref_count;
  my $ppv = ($query_count-($hit_count_on_query-$hit_count_on_ref)==0)? "big" : $hit_count_on_ref/($query_count-($hit_count_on_query-$hit_count_on_ref));
  my $bp_sen = $bp_hit/$ref_bp_count;
  my $bp_ppv = ($query_bp_hit !=0)? $bp_hit/$query_bp_hit : "big";

 print "\tsen=$sen, ppv=$ppv, ref_count=$ref_count, hit_on_ref=$hit_count_on_ref, hit_count_on_query=$hit_count_on_query, query=$query_count\n";
 print "\tbp_sen=$bp_sen, bp_ppv=$bp_ppv, $ref_bp_count, hit=$bp_hit, query=$query_bp_hit\n";
}


