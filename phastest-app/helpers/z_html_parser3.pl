#!/usr/bin/perl -w

package MyParser1;
    	use base qw(HTML::Parser);
    	use LWP::Simple ();
	my $read_flag =0;
	my $sp1 ='';
	my $name ='';
	my $name_flag = 0;
	my $td_flag = 0;
	my $tr_flag = 0;
	my $plasmid_name = '';
	my $acc_num = '';
	my $plasmid_flag =0;
	my $acc_flag =0;
	my $gen_acc_flag =0;
	my $gen_acc_num='';
	
	sub start {
		my ($self, $tagname, $attr, $attrseq, $text) = @_;
		if ($tagname eq 'table' ){
			my $class = $attr->{ class };
			if ($class eq "jig-ncbigrid"){
				$read_flag = 1;
			}
		}elsif($tagname eq 'tr' && $read_flag ==1){
			$plasmid_name = '';
			$acc_num = '';
			$plasmid_flag =0;
			$acc_flag =0;
			$td_flag = 0;
			$tr_flag = 1;
			$gen_acc_num='';
			$gen_acc_flag =0;
		}elsif($read_flag==1 && $tagname eq 'td'){
			$td_flag = 1;			
			if ($plasmid_flag ==1){
				$acc_flag = 1;
				$plasmid_flag = 0;
			}elsif ($acc_flag ==1){
				$gen_acc_flag =1;
				$acc_flag = 0;
			}
		}
		
	}
	sub text{
		my ($self,$text) = @_;
		if ($text eq "Name:"){
			$name_flag =1;
		}elsif ($name_flag ==1){
			$sp1 = $text;
			$name_flag =0;
		}
		if ($td_flag ==1 &&  $text=~/Plasmid/i) {
			$plasmid_name=$text;
			$plasmid_flag =1;
		}elsif($td_flag ==1 && $acc_flag ==1){
			$acc_num = $text;
		}elsif($td_flag ==1 && $gen_acc_flag ==1){
			$gen_acc_num = $text;
		}	

	}
	sub end{
		my ($self, $tagname)=@_;
		if ($tagname eq "table" && $read_flag ==1){
			$read_flag =0;
		}
		elsif($read_flag==1 && $tagname eq 'td'){
			$td_flag =0;
		}
		elsif($tr_flag==1 && $tagname eq 'tr'){
			$tr_flag =0;
			if ($gen_acc_flag ==1 or $acc_flag ==1 ){
				$gen_acc_flag = 0;
				$acc_flag =0;
				if ($gen_acc_num eq ''){
					$gen_acc_num = $acc_num;		
				}
				if ($acc_num eq ''){
					$acc_num= $gen_acc_num;
				}
				print  "-\t-\t$sp1, $plasmid_name\t-\t-\t-\t-\t-\t-\t$gen_acc_num\t$acc_num\n"; 
			
			} 
		}
	}


package MyParser;
    use base qw(HTML::Parser);
    use LWP::Simple ();
    
    my $found_RefSeq_PID= 0;
    my $found_tr = 0;
    my $found_td = 0;
    my $plasmid_num = '';
    my $link = '';    
    my %already_parsed;
    my $t;
    my $td_count=0;
    my $link_flag = 0;
    my $parser1 = MyParser1->new;
    $parser1->handler(start =>  "start", 'self, tagname, attr, attrseq, text' );
	 $parser1->report_tags("table", "tr", "td", "span");
    sub start {
	my ($self, $tagname, $attr, $attrseq, $text) = @_;
	if ($tagname eq 'tr' && $found_RefSeq_PID==1) {
	    $found_tr=1;
	    $td_count=0;
	    $plasmid_num = '';
	    $link = '';
	    $link_flag = 0;
	   
	}elsif ($tagname eq 'td' && $found_tr==1){
		$found_td = 1;
		$td_count++;
	}elsif ($found_tr==1 && $found_td ==1 && $td_count==1 && $tagname eq 'a'){
		$link = $attr->{ href };
		if ($link eq 'http://www.ncbi.nlm.nih.gov/bioproject/0'){
			$link = '';
			$link_flag = 0;
		}else{
			$link_flag = 1;
		}
	}elsif ($found_tr==1 && $found_td ==1 && $td_count==3 && $tagname eq 'a' && $link eq ''){
		$link = $attr->{ href };
		$link_flag = 1;
	}
	
    }
    sub text
       {
         my ($self,$text) = @_;
	if ($text eq 'RefSeq PID'){
		$found_RefSeq_PID= 1;
	}elsif( $found_tr==1 &&  $found_td ==1){
		if ($text eq '' or $text eq '&nbsp;'){
			$text='';
		}
		$text=~s/^[\n\t\s]*//g;
		$text=~s/[\n\t\s]*$//g;
		
		if ($text ne ''){
			$t=$text;
		}
		
	}
	
       }
    sub end{
	my ($self, $tagname)=@_;
	if ($tagname eq 'table' && $found_RefSeq_PID==1){
		$found_RefSeq_PID= 0;
	}elsif($link_flag ==1 && $tagname eq 'a'){
		$link_flag = 0;
	}elsif ($tagname eq 'tr' && $found_tr==1){
		print "\n";
		$found_tr=0;
		#if ($plasmid_num =~/\d+/){
			my $html1 = '';
			$html1 = LWP::Simple::get($link);
    			$parser1->parse( $html1 );
		#}
	}elsif ($found_td ==1 && $tagname eq 'td'){
		$found_td = 0;
		if ($td_count==9){
			$plasmid_num = $t;
		}
		$t = '-' if($t eq '');
		print "$t\t";
		$t='';
	}

   }

package main;
	use LWP::Simple ();
    my $start_time = time;
    my $html = LWP::Simple::get("http://www.ncbi.nlm.nih.gov/genomes/lproks.cgi");
    my $parser = MyParser->new;
     $parser->handler(start =>  "start", 'self, tagname, attr, attrseq, text' );
	 $parser->report_tags("table", "tr", "td", "a");
    $parser->parse( $html );
    my $end_time = time;
    print STDERR  "parse http://www.ncbi.nlm.nih.gov/genomes/lproks.cgi run time: ".($end_time - $start_time)."\n";

