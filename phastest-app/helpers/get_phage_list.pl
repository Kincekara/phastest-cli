#!/usr/bin/perl -w

our @array=();
our $BASE_ROOT = "/var/www/html/phast/current/public";

package MyParser;
    use base qw(HTML::Parser);
    use LWP::Simple ();

    sub start {
	my ($self, $tagname, $attr, $attrseq, $text) = @_;
       }
    sub text
       {
         my ($self,$text) = @_;
	if ($text =~/NC_/){
		push @array , $text;
	}
	
       }
    sub end{
	my ($self, $tagname)=@_;
   }

package main;
	use LWP::Simple ();
	my $tmp_dir ="$BASE_ROOT/tmp_virus_list" ;
	if (!(-d $tmp_dir)){
		mkdir $tmp_dir;
	}
	my $exec_dir = "/apps/phast/project/cgi-bin";
	my $list_filename ="$tmp_dir/phage_list";
    my $html = LWP::Simple::get("http://www.ncbi.nlm.nih.gov/genomes/GenomesGroup.cgi?taxid=10239&opt=Virus&sort=genome");

    my $parser = MyParser->new;
     $parser->handler(start =>  "start", 'self, tagname, attr, attrseq, text' );
    $parser->parse( $html );

	print "From NCBI, total number of virus genome = ".(scalar @array)."\n";
    #check if different with the old phage list file
	my @add_array =();
	if (!(-e $list_filename)){
		system("touch $list_filename");
	}
	my $phage_list_file_content=`cat $list_filename`;
  	foreach my $NC (@array){
		if ($phage_list_file_content !~ /$NC/s){
			push @add_array , $NC; 
		}
	}
	if (scalar @add_array ==0){
		print "There is no extra new virus entry from NCBI\n\n";
		exit;
	}else{
		print "There are ".(scalar @add_array)." new virus entries from NCBI.\n\n";
	}
    open (OUT, ">>$tmp_dir/vgenome.tbl") or die "Cannot append $tmp_dir/vgenome.tbl";
    my $c =0;
    foreach my $NC (@add_array){
	$c++;
	print "Working on $NC, count $c\n";
	system("echo $NC >> $list_filename");
	system("perl $exec_dir/get_gbk.pl  $NC  $tmp_dir");# create .gbk file
          system("perl $exec_dir/gbk2faa.pl $tmp_dir/$NC.gbk"); #create .faa filed
	my $data = `cat $tmp_dir/$NC.gbk`;
	my $legth = '';
	my $sp = '';
	if ($data =~/LOCUS\s+\S+\s+(\d+)\s+bp.*?\nSOURCE\s+(.*?)\n/s){
		$length=$1;
		$sp = $2;
		$sp =~s/\(.*\)//;
	}
	my $str ='';
	my $count=0;
	open (IN, "$tmp_dir/$NC.faa") or die "Cannot open $NC.faa";
	while (my $l = <IN>){
		if ($l=~/>gi\|(\d+)\|/){
			$str .= "\t$1";
			$count++;
		}
	}
	close IN;
	if ($count !=0){	
		print OUT "$NC\t$sp\t$length\t$count$str\n";
	}
	unlink "$tmp_dir/$NC.gbk", "$tmp_dir/$NC.faa";
    }
    close OUT;
	system("touch $tmp_dir/vgenome.tbl; cp -f $tmp_dir/vgenome.tbl $BASE_ROOT/database/vgenome.tbl");
	print "Cp vgenome.tbl to $BASE_ROOT/database/\n";

    exit;


