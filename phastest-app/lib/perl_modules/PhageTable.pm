#!/usr/bin/perl
use strict;
package PhageTable;
use Gene;
use Prophage;
use Record;
use Dictionary;
use IS;

my $IS = new IS();
sub new
{
	my $class = shift;
	my %GI_acc = ();
	my %acc_name = ();
	my %acc_size = ();
	my %acc_num_gene = ();
	my %localacchash = ();
	my %acc_GIs=();
	my %acc_PHAGE=();
	my $self = {
		# constant data
		_file => shift,		# phage genome file
		_GI_acc => \%GI_acc,		# hash GI -> its host's acc
		_acc_name => \%acc_name,	# acc -> genome name
		_acc_size => \%acc_size,		# acc -> genome size
		_acc_num_gene => \%acc_num_gene,	# acc -> number of genes in genome
		_acc_GIs=> \%acc_GIs,
		_acc_PHAGE => \%acc_PHAGE,
		# assessment data
		_localacchash => \%localacchash,	# store informations of a given phage
		_integrase => 0,
        _integrase_names=>'',
        _transposase => 0,
        _transposase_names=>'',
        _structural => 0,
        _structural_names=>'',
		_tags=>'',
		_totalgene => 0,  # for completeness score
		_totalbase => 0,  # for completeness score
		_useIS => shift,	# wether to use IS filter
		_IS => 0			# IS element to exclude
	};
	bless $self, $class;
	
	# load file into constant hashes
	open (R, $self->{_file}) or die "cannot read phage genome file",  $self->{_file};
	while (my $line=<R>){
		chomp($line);# added by Jack on Jan 24, 2013; no this line, the last element when arrayed will be with '\n';
		my @tokens = split(/\t/, $line);
		if ($#tokens > 3){
			my ($acc, $name, $size, $ngene);
			$acc = $tokens[0];
			$name = $tokens[1];
			$size = $tokens[2];
			$ngene = $tokens[3];
		 	next if (defined $acc_name{$acc});	
			$acc_name{$acc} = $name;
			$acc_size{$acc} = $size;
			$acc_num_gene{$acc} = $ngene;
			for my $i (4 .. $#tokens){

				if (defined $GI_acc{$tokens[$i]}){
					#warn "Warning: GI:tokens[$i] has been found in more than one place in phage table";
				}
				else{
					$GI_acc{$tokens[$i]} = $acc;
				}
				$acc_GIs{$acc} .= "$tokens[$i] ";
			}
		}
	}
	return $self;
}

# reset the assessment hash
sub clear {
	my ($self) = @_;
	my %hash = ();
	$self->{_localacchash} = \%hash;
	$self->{_integrase} = 0;
    $self->{_integrase_names}='';
    $self->{_transposase} = 0;
    $self->{_transposase_names}='';
    $self->{_structural} = 0;
    $self->{_structural_names}='';
	$self->{_pre_annotated_phages} = 0;
	$self->{_tags}='';
	$self->{_totalgene} = 0;
	$self->{_totalbase} = 0;
	$self->{_IS} = 0;
}

sub printParts {
	my ($self) = @_;
	my $gsize = 1000000;

	my $int_score = 0;
	my $tra_score = 0;
	my $str_score = 0;
	my $phg_score = min($self->{_pre_annotated_phages} * ($gsize / 5), $gsize);
	my $CDS_score = 0;

	$int_score = min($self->{_integrase} * ($gsize / 20), ($gsize / 4)) if ($self->{_integrase} >= 5);
	$tra_score = min($self->{_transposase} * ($gsize / 60), ($gsize / 5)) if ($self->{_transposase} >= 10);
	$str_score = min($self->{_structural} * ($gsize / 15), $gsize) if ($self->{_structural} >= 5);
	$CDS_score = min(250000 + ($self->{_tags} - 15) * ($gsize / 20), ($gsize * 3 / 4)) if $self->{_tags} >= 15;

	print "Integrase:   ", $self->{_integrase},"\n";
	print "Transposase: ", $self->{_transposase},"\n";
	print "Other keys:  ", $self->{_structural},"\n";
	print "Phage keys:  ", $self->{_pre_annotated_phages},"\n";
	print "IS        :  ", $self->{_IS},"\n";
	print "CDS       :  ", $self->{_tags}, "\n";
	print "Integrase score:   ", $int_score,"\n";
	print "Transposase score: ", $tra_score,"\n";
	print "Other keys score:  ", $str_score,"\n";
	print "Phage keys score:  ", $phg_score,"\n";
	print "CDS score:         ", $CDS_score,"\n";
}

# This subroutine basically adds the prophage region $prophage to the prophage region
# already represented by this PhageTable object, thereby "growing" the prophage region.
# This subroutine then updates the various stats used for scoring the new combined
# prophage region.
# Note: $n and $c are only used for debug statements.
sub assessprophage {
	my ($self, $prophage, $n, $c) = @_;

	$self->{_tags}.=$prophage->getNumberOfCDS();
	$self->{_totalbase} += $prophage->getTail()->getEnd() - $prophage->getHead()->getStart() + 1;
	my $gene_c=0;
	for (my $gene = $prophage->getHead(); $gene ne ''; $gene = $gene->next()){
		$gene_c++;
		if ($self->{_useIS} > 0 and $IS->isExc($gene->getLocalGI()) == 1){
			$self->{_IS} = -1;
		}
		elsif ($self->{_useIS} > 0 and $IS->isInc($gene->getLocalGI()) == 1){
			$self->{_IS} = 1;
		}	
#		$gene->printInfo(0);
#		if ($gene->getProduct() =~m/integrase/i or $gene->getProduct() =~m/site[\s\-]specific recombinase/i or $gene->getProduct() =~m/phage recombinase/i){
		if ($gene->getProduct() =~m/integrase|specific recombinase|phage recombinase/i){
			$self->{_integrase} += 1;
			$self->{_integrase_names}.=$gene->getProduct()."gi".$gene->getLocalGI.", ";	
#			print $gene->getProduct() . "	; gi: ".$gene->getLocalGI."\n";
		}elsif ($gene->getProduct() =~m/transposase/i){
			$self->{_transposase} += 1;
			$self->{_transposase_names}.=$gene->getProduct()."gi".$gene->getLocalGI.", ";
#			print $gene->getProduct() . "	; gi: ".$gene->getLocalGI."\n";
		}elsif ($gene->getProduct() =~m/(capsid|holin|head|coat|envelope|virion|flippase|host|injection|tail|plate|portal|neck|scaffold|structural|antirespressor|anti-repressor)/i){
			$self->{_structural} += 1;
			$self->{_structural_names}.=$gene->getProduct().$gene->getLocalGI.", ";	
#			print $gene->getProduct() . "	; gi: ".$gene->getLocalGI."\n";
		}
		if ($gene->getProduct() =~m/prophage|bacteriophage/i) {
			$self->{_pre_annotated_phages}++;
		}
		
		my $records = $gene->getBLASTresult();
#		print STDERR "aaaa row $n column $c, rec count=".(scalar @$records)."\n" if ($records ne '' && $records ne 'keyword' && $n >=0);
		my $gi_str='';
		if ($records ne '' and $records ne 'keyword'){
			foreach (@$records){
				#if ($_->getDefinition() =~m/integrase/i or $_->getDefinition() =~m/site.specific recombinase/i or $_->getDefinition() =~m/phage recombinase/i)
				if ($_->getDefinition() =~m/integrase|specific recombinase|phage recombinase/i){
					$self->{_integrase} += 1;
					$self->{_integrase_names}.=$_->getDefinition()."ref".$_->getRefAcc.", ";
#					print STDERR "aaaa row $n column $c, integrase=". $_->getGI."\n" if (($n==1||$n==2)&&$n >=0);
				}
				elsif ($_->getDefinition() =~m/transposase/i){
					$self->{_transposase} += 1;
					$self->{_transposase_names}.=$_->getDefinition()."ref".$_->getRefAcc.", ";
#					print STDERR "aaaa row $n column $c transposase=". $_->getGI."\n" if (($n==1||$n==2)&&$n>=0);
				}
				elsif ($_->getDefinition() =~m/(capsid|holin|head|coat|envelope|virion|flippase|host|injection|tail|plate|portal|neck|scaffold|structural|antirespressor|anti-repressor)/i){
					$self->{_structural} += 1;
					$self->{_structural_names}.=$_->getDefinition()."ref".$_->getRefAcc.", ";	
#					print STDERR "aaaa row $n column $c, _strcutural detected\n" if ($n>=0);
				}else{
#					print STDERR "aaaa row $n column $c,no inte , no trans. no struct detected\n" if ($n>=0);
				}
				
				# Here we get the accession number associated with the GI number of the BLAST hit.
				#
				# _acc_PHAGE:
				#      keys: all accession numbers of all viruses from all BLAST hits in this
				#         putative prophage region.
				#      values: the species name for the given accession number.
				# _localacchash is a 2-level hash:
				#      keys: all accession numbers of all viruses from all BLAST hits in this
				#         putative phage region.
				#      values: hash where keys are all (unique) GI numbers of all BLAST hits in
				#         this prophage region from the virus with the accession number in
				#         question, and the values are the corresponding BLAST e-values.
				#         Note: The e-values are only used for breaking ties when assigning an
				#         accession to a prophage region.
				my $ref = $_->getRefAcc();
				my $evalue = $_->getEvalue();
				my ($acc) = $_->getSpecies() =~m/.*?(\w{2}_[^_]+)$/;
				if (defined $self->{_acc_GIs}->{$acc}){
					# Jack, July 20, 2018. Here is not right. Based on GC to find acc 
					# from table DB/vgenome.tbl, there are more accs have same gis.
#					$acc = $self->{_GI_acc}->{$gi};
					$self->{_acc_PHAGE}->{$acc}=$_->getSpecies();
					if (defined $self->{_localacchash}->{$acc}){
						$self->{_localacchash}->{$acc}->{$ref} = $evalue;
					}
					else{
						my %newhash = ();
						$newhash{$ref} = $evalue;
						$self->{_localacchash}->{$acc} = \%newhash;
					}
				}
				else{
	
				}
				my $query_gi= $_->getRefAcc;
				if ($gi_str !~/$query_gi/){
					$self->{_totalgene}++ ;
					$gi_str.="$query_gi ";
				}
			}
		}
		if ($gene == $prophage->getTail()){
			last;
		}
	}

}

# Notes by David Arndt:
# This subroutine returns the overall score for this putative prophage region.
# It also counts the number of unique BLAST hits (gi numbers) for each viral accession
# number (generally species or strain) in the putative prophage region -- considering all
# BLAST hits of all identified viral genes, not just the top hits -- and returns this data
# in $PHAGE_str.
#
# Note: scan.pl calls this subroutine with $gsize == 1,000,000.
#
# Note: For some odd reason, scan.pl calls this subroutine with $compscore == the maximum
# length of a phage region after joining (huh?) == 150,000. I guess it saves memory from
# not having to store that number in two variables!
sub evaluate {
	my ($self, $gsize, $steppenalty, $compscore) = @_;
	# penalty will be total gaps between joined prophages
	# print STDERR "nnn steppenalty=$steppenalty\n";
	my $score = 0 - abs($steppenalty);

	# base score makes sure the phage has integrase and structural protein (transposase+structural for mu-likes)
	if ($self->{_IS} == 1){
		$score += $gsize;
	}
	elsif ($self->{_IS} == -1){
		$score += 0;
	}
	elsif (	($self->{_integrase} > 0 and $self->{_structural} > 0) or 	# regular phage
			($self->{_transposase} > 0 and $self->{_structural} > 0)){ 	# Mu-like phage
		$score += ($gsize / 2);
	}

	# Complete phage synteny.
	if ((index($self->{_structural_names}, "tail") != -1) and ((index($self->{_structural_names}, "head") != -1) or (index($self->{_structural_names}, "capsid") != -1))) {
		$score += ($gsize / 4);
	}

	# High numbers of integrase/transposase and/or structural proteins.
	$score += min($self->{_integrase} * ($gsize / 20), ($gsize / 4)) if ($self->{_integrase} >= 5);			# Max 250000 at 5 integrases.
	$score += min($self->{_transposase} * ($gsize / 60), ($gsize / 5))  if ($self->{_transposase} >= 10);	# Max 200000 at 15 transposases.
	$score += min($self->{_structural} * ($gsize / 15), $gsize) if ($self->{_structural} >= 5);				# Max 1000000 at 15 structural proteins.

	# Numbers of pre-annotated proteins.
	$score += min($self->{_pre_annotated_phages} * ($gsize / 5), $gsize);	# Max 1000000 at 5 pre-annotated phage proteins.

	# Numbers of CDS annotated in general (includes various phage hypothetical proteins).
	$score += min(250000 + ($self->{_tags} - 15) * ($gsize / 20), ($gsize * 3 / 4)) if $self->{_tags} >= 15;	# Max 750,000 at 25 CDS.

	my $hit_names = $self->{_integrase_names}.$self->{_transposase_names}. $self->{_structural_names};
	# phage genome bonus join two prophages if they can make a complete known phage genome
	my $max = 0; # record the max score in the while loop
	my $NC_max='NA';
	my $percentage_max=-1;
	my $hit_count_max=-1;
	my $origin_length_max=-1;
	my $best_avg_log_evalue=9999;
	my @arr=();

	while (my ($acc, $hash) = each %{$self->{_localacchash}}){
		my $curkey = scalar(keys %$hash); # number of genes found
		# print STDERR "$acc, $self->{_acc_PHAGE}->{$acc}, $curkey\n";
		my $new_acc_phage= ACC_PHAGE->new($acc, $self->{_acc_PHAGE}->{$acc}, $curkey);
		push @arr, $new_acc_phage;
		my $value = join " ", keys %$hash;
		my $orikey = $self->{_acc_num_gene}->{$acc};	# number of genes in phage genome (in our database)
		my $avg_log_evalue = 0;
		my $s = 0;

		# Notes by David Arndt:
		# Here $curkey = number of unique BLAST hits (gi numbers, out of all BLAST hits for
		#                all matched viral genes in this putative prophage region) that are
		#                from the phage with accession $acc.
		#      $orikey = number of genes that the phage with accession $acc has according to
		#                our database.
		# If >= 50% of the genes in a small phage (according the our database) are matched by
		# our BLAST hits, or >= 15% of the genes in a large phage, we calculate a bonus to the
		# phage region score. But we only add the bonus for the phage accession number that
		# yields the highest bonus to the putative prophage region's score.
		my $percentage;
		if ($orikey == 0) {
			$percentage = 0;
		}else{
			$percentage= $curkey/$orikey;
		}

		# Find the average log10(e-value) of all BLAST hits for the accession number.
		# This is only used for breaking ties when determining which accession to assign to
		# a candidate prophage region.
		my $sum_log_evalue = 0;
		my $count_log_evalue = 0;
		for my $evalue (values %$hash) {
			if ($evalue > 0.0) {
				$sum_log_evalue += log($evalue)/log(10);
			} else {
				$sum_log_evalue += -324.0; # Treat as e-value of 1e-324, which is about the minimum. https://www.biostars.org/p/297267/
			}
			$count_log_evalue += 1;
		}
		$avg_log_evalue = $sum_log_evalue/$count_log_evalue;

		# Check several metrics to determine the best accession to assign to the candidate
		# prophage region:
		#   First compare the "percentage" score.
		#   If that ties, check the total number of genes.
		#   If that ties, check the average log10(e-value) of BLAST hits.
		#   If that ties, use the accession that is the lowest alphabetically.
		# The reason for breaking ties is to avoid non-deterministic behaviour that would
		# result in different accessions being assigned to a prophage between different runs,
		# even with identical input. If we were to only check ($percentage >= $percentage_max),
		# this non-deterministic behaviour would show up. The reason is that Perl uses a
		# different hash seed each time it is run, for security reasons, resulting in hash
		# elements being processed in different order between different Perl runs.
		if ($percentage >= $percentage_max){
			my $select_new_max = 0;
					# Whether to select the current accession as the new best-scoring accession in
					# place of the old one.
			if ($percentage == $percentage_max){
				# Tied percentage score.
				if ($curkey > $hit_count_max) {
					# Break tie based on which accession has more total genes.
					$select_new_max = 1;
				} else {
					if ($avg_log_evalue < $best_avg_log_evalue) {
						# Break tie based on which accession got a lower average log10(e-value) across
						# all its BLAST hits in the candidate prophage region.
						$select_new_max = 1;
					}
					elsif (($acc cmp $NC_max) == -1) {
						# Break tie by choosing accession that is first alphabetically (generally,
						# therefore, with the lowest accession number).
						$select_new_max = 1;
					}
				}
			}
			else { # $percentage > $percentage_max)
				$select_new_max = 1;
			}

			if ($select_new_max) {
				$percentage_max=$percentage;
				$hit_count_max= $curkey;
				$origin_length_max= $orikey;
				$NC_max=$acc;
				$best_avg_log_evalue=$avg_log_evalue;
			}
		}

		if ($orikey < 20){
			if ($percentage >= 0.5){
				$s = $compscore*$percentage;
			}
		}
		else{
			if ($percentage >= 0.15){
				$s = $compscore*$percentage;
			}
		}
		# print STDERR "pppp $acc, $curkey,'$value', $orikey, $compscore, $s\n";
		if ($s > $max){
			$max = $s;
		}
		# print STDERR "$acc  $curkey / $orikey: $s\n";
	}
	$score += $max;
	@arr= sort {$b->{_count} <=> $a->{_count}} @arr;
	my $PHAGE_str='';
	foreach my $o (@arr){
		$PHAGE_str.=$o->{_PHAGE}."(".$o->{_count}."),";
	}
	$PHAGE_str=~s/,$//s;
	# print STDERR "final = $score, percentage_max=$percentage_max, hit_count_max=$hit_count_max, origin_length_max=$origin_length_max, NC_max=$NC_max\n";
	return $score, $percentage_max, $hit_count_max, $origin_length_max, $NC_max, $PHAGE_str;
}

# assign a completeness (was called viablity) score to a prophage (wrapped in PhageTable class)
# the output is a string (so Jack doesn't need to change his parser)
# output (according to paper): case A: =150=PHAGE NAME=
# case B: =some score=PHAGE NAME=
# case C: =-1=
# the script doesn't calculate score for case C (it will be added by Jack's script) 
sub completeness {
	my ($self) = @_;
	my $mostcompletephage = ''; # the phage that has most genes in the prophage region (by %)
	my $mostpercentage = 0; # percentage of above
	my $most_curkey='';
	my $most_orikey='';
	while (my ($acc, $hash) = each %{$self->{_localacchash}}){
		my $curkey = scalar(keys %$hash);				# number of genes found
		my $orikey = $self->{_acc_num_gene}->{$acc};	# number of genes in phage genome
		if ($curkey == $orikey){ # case A
			return "=150=", $self->{_acc_name}->{$acc},"=";
		}
					my $p;
		if ($orikey ==0){
			$p = 0;
		}else{
			$p = $curkey/$orikey;
		}
		if ($p > $mostpercentage){
			$mostpercentage = $curkey/$orikey;
			$mostcompletephage = $acc;
			$most_curkey=$curkey;
			$most_orikey=$orikey;
		}
	}
	if ($mostpercentage >= 0.5){
	#	my $score = 100*($self->{_totalbase})/($self->{_acc_size}->{$mostcompletephage});
		my $score = int(150*$mostpercentage);
		if ($score > 150){
			$score = 150;
		}
		return "=$score=".($self->{_acc_name}->{$mostcompletephage})."=";
	}
	return "=-1="; # default result for case C
}

sub min {
	my ($x, $y) = @_;
	return $y if ($x >= $y) ;
	return $x;
}

package ACC_PHAGE;
sub new {
	my $class=shift;
	my $self= {
		_acc=>shift,#accession number
		_PHAGE=>shift,#PHAGE name
		_count=>shift #gene count
		};
	bless $self, $class;
	return $self;
};
