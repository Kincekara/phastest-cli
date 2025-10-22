#!/usr/bin/perl
# April 13, 2011 by Rah: handle the sequence file
package DNA;
use strict;
use Prophage;
use Gene;
use Contig;
use Status;
use List::Util qw[min max];

my $MAX_ITER = 50000; # max number of searches (combination)  this prevent very long search (image there'er 100 pairs of matches)
my $cur_iter = 0;
my $integrasedepth = 0.25; # percentatge of genome size the integrase can exist
my $terminalInward = 0.2;
my $minlen = 12;
my $maxlen = 80;
my $terminalOutward = 15000;
my $maxRegionSizeForFullSearch = 100000; # maxium phage size for tamdem phage att search (or it takes forever for exhausive search)

my @contigs = (); # array of Contig objects

sub new
{
	my ($class, $fna) = @_;
    my $self = {
		_definition => "",
		_nucleotide => "",
		_size => 0,
		_globalGC => 0,
		# for contig handling
		_contigsFile => '' # name of contig file. Empty string if no contigs defined.
	};
	open (R, $fna) or die "cannot read fasta file $fna";
	my $seq = "";
	while (my $line =<R>){
		chomp($line);
		if ($line =~m/^>(.+)$/){
			$self->{_definition} = $1;
		}
		elsif ($line =~m/(\w+)/){
			$seq = $seq.$1;
		}
	}
	$self->{_nucleotide} = $seq;
	$self->{_size} = length($seq);
	close (R);
	my $GC = 0;
	while ($seq=~m/[CG]/g){
		$GC++;
	}
	$self->{_globalGC} = $GC;
	bless $self, $class;
    return $self;
}

sub loadContigs {
	my ( $self, $contigsFile ) = @_;
	$self->{_contigsFile} = $contigsFile;
	
	open (R, $contigsFile) or die "cannot open contig positions file $contigsFile";
	my $i = 0;
	while (my $line = <R>){
		chomp($line);
		if ($line =~m/^(.*)\t(\d+)\t(\d+)\t(\d+)$/){
			my $contig = new Contig($2, $3, $1);
			$contig->setId($i);
			push(@contigs, $contig);
			$i++;
		}
	}
	close (R);
}

sub getContigById {
	my ( $self, $id ) = @_;
	
	for (my $i = 0; $i < scalar(@contigs); $i++) {
		my $contig = $contigs[$i];
		if ($contig->getId() == $id) {
			return $contig;
		}
	}
	
	return undef;
}

# Returns the id of the contig at the given position in the overall sequence of concatenated contigs.
sub getContigIdForPos {
	my ( $self, $pos ) = @_;
	
	if ($self->{_contigsFile} eq "") {
		return -1; # this sequence has no defined contigs
	}
	
	for (my $i = 0; $i < scalar(@contigs); $i++) {
		my $contig = $contigs[$i];
		if ($contig->isInContig($pos)) {
			return $contig->getId();
		}
	}
	
	print STDERR "ERROR: Could not find contig at $pos!\n";
	exit(1);
}

sub getSize {
	my ($self) = @_;
	return $self->{_size};
}

sub getGlobalGC {
	my ($self) = @_;
	return $self->{_globalGC};
}

sub getDefinition {
	my ($self) = @_;
	return $self->{_definition};
}

# input starts at 1
sub getGC {
	my ($self, $from, $to) = @_;
	my $str = substr($self->{_nucleotide}, $from-1, $to-$from+1);
	my $GC = 0;
	while ($str=~m/[CG]/g){
		$GC++;
	}
	return $GC*100.0/($to-$from+1);
}
sub constructMatch {
    my $self = {
		_str => shift,
		_pA => shift,
		_pB => shift,
		_len => shift,
		_attL => '',
		_attR => '',
		_score => 0,	# a score used to choose best att
		_inRNA => 0,
		_inORF => 0,
		_nearInt => 0,
		_fit => 0 # absolute difference between cover region and original region
	};
	if ($self->{_pA} < $self->{_pB}){
		$self->{_attL} = $self->{_pA};
		$self->{_attR} = $self->{_pB};
	}
	else{
		$self->{_attL} = $self->{_pB};
		$self->{_attR} = $self->{_pA};
	}
	return $self;
}

sub checkAtt {
	my ($att, $prophage) = @_;
	my $inORF = 0;
	my $inRNA = 0;
	for (my $cursor = $prophage->getHead(); $cursor ne ''; $cursor = $cursor->previous()){
		if ($cursor->getEnd() < $att->{_attL}){
			last;
		}
		if ($cursor->isInGene($att->{_attL}) == 1){
			if ($cursor->isRNA() > 0){
				$inRNA = 1;
			}
			else{
				$inORF = 1;
			}
		}
		if ($cursor->isInGene($att->{_attR}) == 1){
			if ($cursor->isRNA() > 0){
				$inRNA = 1;
			}
			else{
				$inORF = 1;
			}
		}
	}
	for (my $cursor = $prophage->getHead(); $cursor ne ''; $cursor = $cursor->next()){
		if ($cursor->getStart() > $att->{_attR}){
			last;
		}
		if ($cursor->isInGene($att->{_attL}) == 1){
			if ($cursor->isRNA() > 0){
				$inRNA = 1;
			}
			else{
				$inORF = 1;
			}
		}
		if ($cursor->isInGene($att->{_attR}) == 1){
			if ($cursor->isRNA() > 0){
				$inRNA = 1;
			}
			else{
				$inORF = 1;
			}
		}
	}
	$att->{_inORF} = $inORF;
	$att->{_inRNA} = $inRNA;
	$att->{_fit} = abs($prophage->getHead()->getStart() - $att->{_attL})+abs($prophage->getTail()->getEnd() - $att->{_attR});
}

# find att based on phage conditions
# there're three cases:
# 1. (normal phage) single integrase
# 2. (mu like) no integrase but transposase
# 3. (possible piggyback) multiple integrase (may or may not exist transposase)
#
# case 1 requires to find a single pair of attL and attR
# the search will be in the following priority:
# a) if a tRNA is near a terminal (and not exclude the integrase and reads 
# from the terminal to the rest of the phage), search the tRNA against
# the rest of the other end (from 20% of the phage size inward to 15000 bps outward)
# b) if an integrase is near a terminal, search the 400 bps against
# the other end (from 20% of the phage size inward to 15000 bps outward)
# c) neither tRNA or integrase found to be near terminal then the two ends
# (15000 - 20% into the phage) are used to find att
# if more than one att found, the pair that fits the phage the best 
# and do not over lap important gene will be chosen
# 
# case 2 is like c) in case 1, +15000 - 20% into the phage region will be searched for att pairs
# putative att pair must bracket the transposase and the best one best fits the phage
# updated April 19, by Rah: this prediction is removed since the repeat for Mu transposase is usually too smaller
# 
# and 3 require the entire phage +15000 bps of each end to be searched for all
# possible att pairs. putative att pairs must not insert into proteins (at least phage proteins)
# and include at least one integrase or transposase. If mutiple att pairs corresponds to the same
# integrase/transposase then the ones that are within 400 bps of their integrase and/or with tRNA
# then one that contains the most ORF will be chosen. The results may contains
# mutiple att pairs

sub findAttachmentSite {
	my ($self, $prophage) = @_;
	my $result;
# 	print STDERR "search prophage $prophage\n";
	my $int = $prophage->getIntegrase(); # array of integrases (5->3 order)
	my $tran = $prophage->getTransposase(); # array of transposase (5->3 order)
	my $RNA = $prophage->getRNA();# array of RNA (5->3 order)
	my $inputSize = $prophage->getTail()->getEnd() - $prophage->getHead()->getStart();
#	$prophage->printInfo(1);
	my @attachsites = ();
	
#	print STDERR "FIND ATT: ", $prophage->getHead()->getStart(), " - ", $prophage->getTail()->getStart()," int ", $#$int+1, " trans ", $#$tran+1,"\n";
	
	my $contig = undef;
	if ($self->{_contigsFile} ne "") {
		# sequence consists of concatenated contigs
		$contig = $self->getContigById($prophage->getHead()->getContigId());
	}
	
	if ($#$int+1 == 1 or ($#$int+1 > 1 and $inputSize > $maxRegionSizeForFullSearch)){
#	print STDERR "CASE 1: ", $prophage->getHead()->getStart()," ... ", $$int[0]->getStart(), ", ", $$int[0]->getEnd()," ... ",	$prophage->getTail()->getEnd(),"\n"; 
		my $integrase = $$int[0];
		my @rawAtt = ();
		# case 1
		# check whether the integrase is near the terminal of the phage
#		print "INTEGRAE\n";
		
		if (($integrase->getStart() - $prophage->getHead()->getStart())/$inputSize < $integrasedepth){
			# -int-------
# 			print STDERR "int found near 5 end\n";
			$result = $self->search($integrase->getStart()-400,
						$integrase->getStart()-1,
						int($prophage->getTail()->getEnd() - $inputSize*$terminalInward),
						$prophage->getTail()->getEnd()+$terminalOutward,
						$contig);
		}
		elsif (($prophage->getTail()->getEnd() - $integrase->getEnd())/$inputSize < $integrasedepth){
			# -------int-
#			print STDERR "int found near 3 end\n";
			$result = $self->search($integrase->getStart()-400,
						$integrase->getStart()-1,
						$prophage->getHead()->getStart() - $terminalOutward,
						int($prophage->getHead()->getStart() + $inputSize*$terminalInward),
						$contig);
		}
		foreach (@$result){
			&checkAtt($_, $prophage);
			$_->{_nearInt} = 1;
			push @rawAtt, $_;
		}
#		print "tRNA\n";
		foreach (@$RNA){
			my $end;
			if ($_->is5to3() == 1 and $_->getEnd() < $integrase->getStart()){
				# -->RNA>---int---
				$end = int($prophage->getTail()->getEnd() - $inputSize*$terminalInward);
				if ($end < $_->getEnd()){	# the region will be too small
					next;
				}
				$result = $self->search($_->getStart(),
							$_->getEnd()-$minlen,
							$end,
							$prophage->getTail()->getEnd()+$terminalOutward,
							$contig);
			}
			elsif ($_->is5to3() == 0 and $_->getStart() > $integrase->getEnd()){
				# --int---<RNA<---
				$end = int($prophage->getHead()->getStart() + $inputSize*$terminalInward);
				if ($end > $_->getStart()){	# the region will be too small
					next;
				}
				$result = $self->search($_->getStart(),
							$_->getEnd()-$minlen,
							$prophage->getHead()->getStart() - $terminalOutward,
							$end,
							$contig);
			}
		}
		foreach (@$result){
			&checkAtt($_, $prophage);
			push @rawAtt, $_;
		}
#		print "NONE MARKER\n";
		if ($#rawAtt == -1){
			# -----int-----
			my $leftEnd = $integrase->getStart() -1;
			if ($prophage->getHead()->getStart() + $inputSize*$terminalInward < $leftEnd){
				$leftEnd = int($prophage->getHead()->getStart() + $inputSize*$terminalInward);
			}
			my $rightEnd = $integrase->getEnd()+1;
			if ($prophage->getTail()->getEnd() - $inputSize*$terminalInward > $rightEnd){
				$rightEnd = int($prophage->getTail()->getEnd() - $inputSize*$terminalInward);
			}
			$result = $self->search($prophage->getHead()->getStart()-$terminalOutward,
						$leftEnd,
						$rightEnd,
						$prophage->getTail()->getEnd()+$terminalOutward,
						$contig);
		}
		foreach (@$result){
			&checkAtt($_, $prophage);
			push @rawAtt, $_;
		}
		
#		print "Selected the best att\n";
		# select the best att
		my $best = '';
		foreach (@rawAtt){
			if ($_->{_nearInt} > 0){
				$_->{_score} += 2000000;
			}
			if ($_->{_inRNA} > 0){
				$_->{_score} += 1000000;
			}
			if ($_->{_inORF} > 0){
				$_->{_score} -= 100000;
			}
#			print" ATT L ", $_->{_attL}, " ATT R ", $_->{_attR}," score = ", $_->{_score},"\n";
			#$_->{_score} -= $_->{_fit};
			$_->{_score} += length ($_->{_str});
			if ($best eq '' or $best->{_score} < $_->{_score}){
				$best = $_;
	#			print "best socre = ", $best->{_score}, "\n";
			}
		}
		if ($best ne ''){
			push @attachsites, $best;
		}
#		print "CASE1 finished\n";
		# end of case 1
	}
	elsif ($#$int+1 == 0){
# 	print "CASE 2: no att site search will be attempted\n";
		# case 2
		# ---tran---
#		my $result = $self->search($prophage->getHead()->getStart()-$terminalOutward,
#					$$tran[0]->getStart()-$minlen,
#					$$tran[$#$tran]->getEnd()+1,
#					$prophage->getTail()->getEnd()+$terminalOutward);
#		
		# choose the best att
#		my @rawAtt = ();
#		foreach (@$result){
#			&checkAtt($_, $prophage);
#			push @rawAtt, $_;
#		}
#		my $best = '';
#		foreach (@rawAtt){
#			if ($_->{_inORF} > 0){
#				$_->{_score} -= 100000;
#			}
#			$_->{_score} -= $_->{_fit};
#			if ($best eq '' or $best->{score} < $_->{_score}){
#				$best = $_;
#			}
#		}
#		if ($best ne ''){
#			push @attachsites, $best;
#		}
		# end of case 2
	}
	else{
#	print STDERR "CASE 3: \n";
		# case 3 (piggyback multiple integrase multiple att site)
		# search entire phage + 15000 bps pass each terminal
		my $result = $self->fullsearch($prophage->getHead()->getStart()-$terminalOutward,
		$prophage->getTail()->getEnd()+$terminalOutward, $contig);
		
		# pick up the best atts
#		print "we got ",$#$result+1," pairs\n";
		my @filtered = ();
		foreach (@$result){
			&checkAtt($_, $prophage);
			my $valid = 0; # must contains at least one integrase or transposase
			for my $n (0 .. $#$int){
				if ($_->{_attL} + $minlen < $$int[$n]->getStart() and $_->{_attR} > $$int[$n]->getEnd()){
					$valid = 1;
					last;
				}
			}
			if ($valid == 0){
				for my $n (0 .. $#$tran){
					if ($_->{_attL} + $minlen < $$tran[$n]->getStart() and $_->{_attR} > $$tran[$n]->getEnd()){
						$valid = 1;
						last;
					}
				}
			}
			if ($valid == 1 and $_->{_inORF} == 0){ # ORF not include tRNA
#				print "$_->{_str}   $_->{_attL}  $_->{_attR} RNA ($_->{_inRNA}) ORF ($_->{_inORF})\n";
				push @filtered, $_;
			}
		}
		
#		print "check combinations of array of total ", $#filtered+1,"\n";
		# now check all possible combinations of attachment sites
		# this is tracked use an array of binary numbers so every combination
		# can be treated as a string of "01001" (use the 2nd and 5th atts).
		# combinations are evaluate to check whether it is valid
		# att must fully included or exclude by each other and every att should have
		# a corresponding int
		
		# initiate the binary tracker
		my @tracker = ();
		for my $i (0 .. $#filtered){
			$tracker[$i] = 0;
		}
		my @final = ();
		my @combineIntTran = ();
		push @combineIntTran, @$int;
#		push @combineIntTran, @$tran;
#		print "FILTERED = ",$#filtered+1,"\n";
		$cur_iter = 0;
		&simCombination(\@tracker, 0, \@filtered, new Status($prophage, \@combineIntTran), \@final);
#		print "ITER = $cur_iter++\n";
#		print "SIZE OF FINAL = ", $#final+1, "\n";
		# pick up the best attachment site combinations
		my $best = '';
		foreach (@final){
			if ($best eq '' or $best->{_score} < $_->{_score}){
				$best = $_;
			}
		}
		
		if ($best ne ''){
			push @attachsites, @{$best->{_att}};
#			my $aa = $best->{_att};
#			print "best = $best aa = $aa and $#$aa\n";
#			print "SIZE before = ", $#$aa+1, " SIZE after = ", $#attachsites+1,"\n";
		}
#		print "FINAL ", \@final,"\n";
#		exit;
	}
	$prophage->setAttachmentSite(\@attachsites);
#	print STDERR "exist att finding funct\n";
#	exit;
}

sub simCombination {
	my ($tracker, $cur, $atts, $status, $final) = @_;
	if ($cur > $#$tracker or $cur_iter > $MAX_ITER){	# terminal case
#		$status->printInfo();
		$status->{_score} = $status->getScore();
		my $record = new Status("", "");
		$status->cp($record);
		$$final[$#$final+1] = $record;
		if ($cur_iter <= $MAX_ITER){
			$cur_iter++;
		}
		return;
	}
	for my $i (0 .. 1){
		if ($i == 0){
			&simCombination($tracker, $cur+1, $atts, $status, $final);
		}
		elsif ($status->add($$atts[$cur]) == 1){
			$$tracker[$cur] = $i;
			&simCombination($tracker, $cur+1, $atts, $status, $final);
			$status->rm($$atts[$cur]);
		}
	}
	return;
}

# hash-based search for faster speed
# When the sequence consists of contigs concatenated together, $contig will be the Contig
# object where the current phage is located. Otherwise it will be undef.
sub search {
	my ($self, $A5, $A3, $B5, $B3, $contig) = @_;
#	print "search att A($A5, $A3)  B($B5, $B3)\n";

	if (!defined($contig)) {
		if ($A5 < 1){
			$A5 = 1;
		}
		if ($A3 > $self->{_size}){
			$A3 = $self->{_size};
		}
		if ($B5 < 1){
			$B5 = 1;
		}
		if ($B3 > $self->{_size}){
			$B3 = $self->{_size};
		}
	} else {
		my $contig_id = $contig->getId();
		
		if ($A5 < 1 || $self->getContigIdForPos($A5) != $contig_id) {
			$A5 = $contig->getStart();
		}
		if ($A3 > $self->{_size} || $self->getContigIdForPos($A3) != $contig_id){
			$A3 = min($contig->getEnd(), $A3+$terminalOutward);
		}
		if ($B5 < 1 || $self->getContigIdForPos($B5) != $contig_id){
			$B5 = $contig->getStart();
		}
		if ($B3 > $self->{_size} || $self->getContigIdForPos($B3) != $contig_id){
			$B3 = $contig->getEnd();
		}
	}
	my %regionB = ();
	# put region B to hash
	for (my $i = $B5; $i <= $B3-$minlen+1; $i++){
		my $key = substr($self->{_nucleotide}, $i-1, $minlen);
		if ($key =~m/[^ATCG]/i){ # will skip any miss sequence
			next;
		}
		if (! defined $regionB{$key}){
			my @array = ($i);
			$regionB{$key} = \@array;
		}
		else{
			push (@{$regionB{$key}}, $i);
		}
	}
	# loop region A to find hits in the hash
	my @matches =();
	for (my $i = $A5; $i <= $A3-$minlen+1; $i++){
		my $key = substr($self->{_nucleotide}, $i-1, $minlen);
		if ($key =~m/[^ATCG]/i){ # will skip any miss sequence
			next;
		}
#		print "look for '$key'\n";
		if (defined $regionB{$key}){
			# identical regions found, try to extend it
#			print "Identical segment found ($key) at A($i) B(${$regionB{$key}}[0])\n";
			my $finalextension = 0;
			for (my $ext = 1; $ext <= $A3-$i-$minlen; $ext++){
				my $nextA = substr($self->{_nucleotide}, $i-2+$minlen+$ext, 1);
				my $ok = 1;
				foreach (@{$regionB{$key}}){
					my $nextB = substr($self->{_nucleotide}, $_-2+$minlen+$ext, 1);
					if ($nextA ne $nextB){
						$ok = 0;
						last;
					}
				}
				if ($ok == 0){
					last;
				}
#				print "extend to $nextA\n";
				$finalextension = $ext;
			}
			my $finalkey = $key;
			if ($finalextension > 0){
				$finalkey = substr($self->{_nucleotide}, $i-1, $minlen+$finalextension);
			}
#			print "FINAL $key\n";
			if ($minlen+$finalextension > $maxlen){
				next;
			}
			foreach (@{$regionB{$key}}){
				push @matches, &constructMatch($finalkey, $i, $_, $minlen+$finalextension);
			}
#			push (@matches, &constructMatch($key, $i, $regionB{$key}));
			if ($finalextension > 0){
				$i += $finalextension;
			}
		}
	}
	return \@matches;
}

# search entire region
# When the sequence consists of contigs concatenated together, $contig will be the Contig
# object where the current phage is located. Otherwise it will be undef.
sub fullsearch {
	my ($self, $from, $to, $contig) = @_;
	
	if (!defined($contig)) {
		if ($from < 1){
			$from = 1;
		}
		if ($to > $self->{_size}){
			$to = $self->{_size};
		}
	} else {
		my $contig_id = $contig->getId();
		
		if ($from < 1 || $self->getContigIdForPos($from) != $contig_id){
			$from = $contig->getStart();
		}
		if ($to > $self->{_size} || $self->getContigIdForPos($to) != $contig_id){
			$to = $contig->getEnd();
		}
	}
	
#	print "exhaustive search on $from - $to\n";
	my %full = (); # hash the entire region to this variable
	for (my $i = $from; $i <= $to-$minlen+1; $i++){
		my $key = substr($self->{_nucleotide}, $i-1, $minlen);
		if ($key =~m/[^ATCG]/i){ # will skip any miss sequence
			next;
		}
		if (! defined $full{$key}){
			my @array = ($i);
			$full{$key} = \@array;
		}
		else{
			push (@{$full{$key}}, $i);
		}
	}
#	print "done hash now iterate the keys\n";
	# now find all attachment site
	my @matches =();
	for (my $base = $from; $base <= $to-$minlen+1; $base++){
		my $key = substr($self->{_nucleotide}, $base-1, $minlen);
		if ($key =~m/[^ATCG]/i){ # will skip any miss sequence
			next;
		}
#		print "iterates $base..\n";
		my $hits = $full{$key}; # should always find at least one hit;
		my $availableTarget = 0;
		foreach (@$hits){
			if ($_ <= $base){ # don't care base itself and ones before base
				next;
			}
			$availableTarget++;
		}
		if ($availableTarget == 0){ # nothing left to pair
			next;
		}
		my $finalextension = 0;
		for (my $ext = 1; $ext <= $to-$base-$minlen and $ext <= $maxlen+1; $ext++){
			my $nextA = substr($self->{_nucleotide}, $base-2+$minlen+$ext, 1);
			my $ok = 1;
			foreach (@$hits){
	#			print "base $base  $$hits[0] and $$hits[1]\n";
				if ($_ <= $base){ # don't care base itself and ones before base
					next;
				}
				my $nextB = substr($self->{_nucleotide}, $_-2+$minlen+$ext, 1);
				if ($nextA ne $nextB){
					$ok = 0;
					last;
				}
			}
			if ($ok == 0){
				last;
			}
# 			print STDERR "extend to $nextA  ($ext)\n";
			$finalextension = $ext;
		}
		my $finalkey = $key;
		if ($finalextension > 0){
			$finalkey = substr($self->{_nucleotide}, $base-1, $minlen+$finalextension);
		}
#		print "FINAL '$finalkey' at $base ext $finalextension to ",$$hits[$#$hits], "\n";
		if ($minlen+$finalextension <= $maxlen){
			foreach (@$hits){
				if ($_ <= $base){ # don't care base itself and ones before base
					next;
				}
				push @matches, &constructMatch($finalkey, $base, $_, $minlen+$finalextension);
			}
		}
		
		if ($finalextension > 0){
			$base += $finalextension;
		}
		
	}
#	print "has found all ", $#matches+1," possible pairs\n";
	return \@matches;
}

1;
