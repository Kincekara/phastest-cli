#!/usr/bin/perl
# April 11, 2011 by Rah. Class to store predicted prophage
use strict;
use Gene;
package Prophage;

sub new
{
	my $class = shift;
    my $self = {
		# constant variables
		_from => "", 		# the first Gene (5->3) 
		_to => "",			# the last Gene (5->3)
		_geneHash => "",	# a hash keyed by Gene object reference
		_integrase => "",
		_att => '',
		_p5 => 0,
		_p3 => 0,
		_GC => 0,
		_completeness => "=-1=", # used by completeness score 
		# linkers
		_next => "",
		_previous => "",
		_percentage_max=>-1, # in digital mumber
                _hit_count_max=>-1,
                _origin_length_max=>-1,
                _NC_max=>'',
		_NC_acc=>'',
		_NC_PHAGE=>''

	};
	my %hash = ();
	$self->{_geneHash} = \%hash;
	bless $self, $class;
    return $self;
}

sub set_percentage_max{
	my ($self, $p)=@_;
	$self->{_percentage_max}=$p;
}
sub get_percentage_max{
	my $self=shift;
	return $self->{_percentage_max};
}
sub set_hit_count_max{
	my ($self, $h)=@_;
	$self->{_hit_count_max}=$h;
}
sub get_hit_count_max{
	my $self=shift;
	return $self->{_hit_count_max};
}
sub set_origin_length_max{
	my ($self, $o)=@_;
	$self->{_origin_length_max}=$o;
}
sub get_origin_length_max{
	my $self=shift;
	return $self->{_origin_length_max};
}
sub set_NC_PHAGE{
	my $self=shift;
	$self->{_NC_PHAGE}=shift;
}
sub get_NC_PHAGE{
	my $self=shift;
	return $self->{_NC_PHAGE};
}
sub set_NC_max{
	my ($self, $NC)=@_;
	$self->{_NC_max}=$NC;
}
sub get_NC_max{
	my $self=shift;
	return $self->{_NC_max};
}
sub set_NC_acc{
	my ($self, $acc)=@_;
	$self->{_NC_acc}=$acc;
}
sub get_NC_acc{
	my $self=shift;
	return $self->{_NC_acc};
}
sub setCompleteness {
	my ($self, $completeness) = @_;
	$self->{_completeness} = $completeness;
}
sub getCompleteness {
	my ($self) = @_;
	return $self->{_completeness};
}
sub getNumberOfCDS {
	my ($self) = @_;
	my $n = 0;
	for (my $cur = $self->getHead(); $cur ne ''; $cur = $cur->next()){
		$n++;
		if ($cur == $self->getTail()){
			last;
		}
	}
	return $n;
}
sub calculateGC {
	my ($self, $sequence) = @_;
	if ($self->{_p5} != 0 and $self->{_p3} != 0){
		$self->{_GC} = $sequence->getGC($self->{_p5}, $self->{_p3});
	}
}
sub get5end {
	my ($self) = @_;
	return $self->{_p5};
}
sub get3end {
	my ($self) = @_;
	return $self->{_p3};
}
sub getSize {
	my ($self) = @_;
	return $self->{_p3} - $self->{_p5};
}
sub getHeadTailSize {
	my ($self) = @_;
	return $self->getTail()->getEnd() - $self->getHead()->getStart();
}
sub setAttachmentSite {
	my ($self, $att) = @_;
	$self->{_att} = $att;
	my $p5 = $self->getHead()->getStart();
	my $p3 = $self->getTail()->getEnd();
	foreach (@$att){
		if ($_->{_attL} < $p5){
			$p5 = $_->{_attL};
		}
		if ($_->{_attR} > $p3){
			$p3 = $_->{_attR};
		}
	}
	$self->{_p5} = $p5;
	$self->{_p3} = $p3;
}

sub getHead {
	my ($self) = @_;
	return $self->{_from};
}

sub getTail {
	my ($self) = @_;
	return $self->{_to};
}

sub next {
	my ($self) = @_;
	return $self->{_next};
}

sub setNext {
	my ($self, $next) = @_;
	$self->{_next} = $next;
}

sub previous {
	my ($self) = @_;
	return $self->{_previous};
}

sub setPrevious {
	my ($self, $previous) = @_;
	$self->{_previous} = $previous;
}

sub insert {
	my ($self, $gene) = @_;
	my $hash = $self->{_geneHash};
	if (defined $hash->{$gene}){
		return 0;
	}
	$hash->{$gene} = 1;
	if ($self->{_from} eq ''){
		$self->{_from} = $gene;
		$self->{_to} = $gene;
	}
	elsif ($gene->getStart() <= $self->{_from}->getStart()){
		$self->{_from} = $gene;
	}
	elsif ($gene->getStart() >= $self->{_to}->getStart()){
		$self->{_to} = $gene;
	}
}

sub isInclude {
	my ($self, $gene) = @_;
	my $hash = $self->{_geneHash};
	if (defined $hash->{$gene}){
		return 1;
	}
	return 0;
}

sub hasIntegrase {
	my ($self) = @_;
	for (my $cursor = $self->getHead(); $cursor->getStart() <= $self->getTail()->getStart(); $cursor = $cursor->next()){
		if ($cursor->isIntegrase() == 1){
			return 1;
		}
		
		# If cursor points to the last element, break out of loop before trying to call getStart() on the empty string (causing a crash).
		if ($cursor->getStart() == $self->getTail()->getStart()) {
			last;
		}
	}
	return 0;
}

sub getIntegrase {
	my ($self) = @_;
	my @int = ();
	for (my $cursor = $self->getHead(); $cursor ne '' and $cursor->getStart() <= $self->getTail()->getStart(); $cursor = $cursor->next()){
		if ($cursor->isIntegrase() == 1){
			push (@int, $cursor);
		}
	}
	return \@int;
}

sub getTransposase {
	my ($self) = @_;
	my @tra = ();
	for (my $cursor = $self->getHead(); $cursor ne '' and $cursor->getStart() <= $self->getTail()->getStart(); $cursor = $cursor->next()){
		if ($cursor->isTransposase() == 1){
			push (@tra, $cursor);
		}
	}
	return \@tra;
}
sub getRNA {
	my ($self) = @_;
	my @rna = ();
	for (my $cursor = $self->getHead(); $cursor ne '' and $cursor->getStart() <= $self->getTail()->getStart(); $cursor = $cursor->next()){
		if ($cursor->isRNA() != 0){
			push (@rna, $cursor);
		}
	}
	return \@rna;
}
sub printInfo {
	my ($self, $flag) = @_;
	print "Prophage region:\n";
	for (my $gene = $self->getHead(); $gene ne ''; $gene = $gene->next()){
		$gene->printInfo($flag);
		if ($gene == $self->getTail()){
			last;
		}
	}
	print "\n";
}
sub structTag {
	my $self = {
		_pos => shift,
		_feature => shift,
		_anno => shift 
	};
	return $self;
}
sub printPhagefinderResult {
	my ($self) = @_;
	print "\n";
	print "END5     FEATNAME           ANNOTATION OF BEST HIT FROM PHAGE DB ANNOTATION, OR PFAM-A HMM HIT\n";
	print "............................................................................................................................................................................\n";
	my $att = $self->{_att};
	my @attall = ();
	if ($#$att == 0){
		push @attall, &structTag($$att[0]->{_attL}, "attL", $$att[0]->{_str});
		push @attall, &structTag($$att[0]->{_attR}, "attR", $$att[0]->{_str});
	}
	else{
		for my $i (1 .. $#$att+1){
			push @attall, &structTag($$att[$i-1]->{_attL}, "attL$i", $$att[$i-1]->{_str});
			push @attall, &structTag($$att[$i-1]->{_attR}, "attR$i", $$att[$i-1]->{_str});
		}
	}
	
	# sort att by index
	for my $x (0 .. $#attall){
		for my $y ($x+1 .. $#attall){
			if ($attall[$y]->{_pos} < $attall[$x]->{_pos}){
				my $temp = $attall[$x];
				$attall[$x] = $attall[$y];
				$attall[$y] = $temp;
			}
		}
	}
	my $it = 0;
	for (my $cur = $self->getHead(); $cur ne ''; $cur = $cur->next()){
		while ($it <= $#attall and $attall[$it]->{_pos} < $cur->getStart()){
			print sprintf("%-11s%-19s%s\n", $attall[$it]->{_pos}, $attall[$it]->{_feature}, $attall[$it]->{_anno});
			$it++;
		}
		my $end = $cur->getStart();
		my $feat = $cur->getLocalGI();
		my $tophit = ''; my $records=$cur->getBLASTresult();
		#print STDERR $self->get_NC_max.":".$self->get_NC_acc."\n";
		if ($records ne '' and $records ne 'keyword'){
			foreach my $r (@$records){
				my $gi=$r->getGI();
				my $gis=$self->get_NC_acc;
				if ($gis=~/$gi/) { # && $self->get_percentage_max()>0.3){
					$tophit=$r;
					last;
				}
			}
			if ($tophit eq ''){
				foreach my $r (@$records){
					my $def = $r->getDefinition();
					if ($def=~m/integrase|specific recombinase|phage recombinase|transposase|capsid|holin|head|coat|envelope|virion|flippase|host|injection|tail|plate|portal/i){
						$tophit=$r;
						last;
					}
				}
			}
			$tophit = $$records[0] if ($tophit eq '');
		}
		# modify to fit the format
		if ($cur->isRNA() == 1){
			$feat = 'tRNA';
		}
		elsif ($cur->isRNA() == 2){
			$feat = 'tmRNA';
		}
		elsif ($cur->isRNA() == 3){
			$feat = 'rRNA';
		}
		if ($cur->is5to3() == 0){
			$end = $cur->getEnd();
		}
		print sprintf ("%-11s%-19s", $end, $feat);
		if ($tophit ne ''){
			unless ($tophit eq "keyword") {
				print sprintf ("gi|%d|ref|%s|, TAG = %s, E-VALUE = %s\n", $tophit->getGI(), $tophit->getRefAcc(), $tophit->getSpecies(), $tophit->getEvalue());
				print sprintf ("%-10s %-10s         ", "", "");
			}
		}
		print sprintf ("[ANNO] %s; %s\n", $cur->getProduct(), $cur->getLtag());
		
		if ($cur == $self->getTail()){
			last;
		}
	}
	while ($it <= $#attall){
		print sprintf("%-11s%-19s%s\n", $attall[$it]->{_pos}, $attall[$it]->{_feature}, $attall[$it]->{_anno});
		$it++;
	}
	
	if ($self->next() eq ''){
		print "............................................................................................................................................................................\n";
	}
	else{
		print "----------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n";
	}
}
1;
