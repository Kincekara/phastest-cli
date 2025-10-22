#!/usr/bin/perl
# April 9, 2011 by Rah
# A class represent a gene (or tRNA) in gbk file
# the data structural of gene object is made to use as double linked list
use strict;
package Gene;

# build a gene (note)
sub new
{
	my $class = shift;
	my $self = {
		# permental data
		_localgi => shift, 	# gi number assigned to the gene
		_start => shift, 	# start index (5->3)
		_end => shift, 		# end index (5->3)
		_strand => shift,	# either '+' 5->3 or '-' 3->5
		_product => shift,	# functional description of the protein or tRNA
		_type => shift,		# either 'p' protein, 't' tRNA, 'tm' tmRNA, or 'r' rRNA
		_ltag => shift,		# lotus tag
		# for contig handling
		_contigId => '', # contig id
		# data for blast
		_namerelated => 0,	# whether the product name is related to phage 
		_blast => 0,		# blast status (0 = never do blast, 1 = will do blast, 2 = done blast)
		_hit => '',			# blast results (array of Record objects, or the string "keyword")
		# DBSCAN data
		_visited => 0,		# flag (visited v.s. unvisited)
		# linkers
		_next => '',
		_previous => ''
	};
	bless $self, $class;
	return $self;
}

sub getLtag {
	my ($self) = @_;
	return $self->{_ltag};
}

sub isVisited {
	my ($self) = @_;
	return $self->{_visited};
}

sub setVisited {
	my ($self, $status) = @_;
	$self->{_visited} = $status;
}

sub setBLASTresult {
	my ($self, $result) = @_;
 	$self->{_hit} = $result;
}

sub getBLASTresult {
	my ($self) = @_;
	return $self->{_hit};
}

sub getBLASTstatus {
	my ($self) = @_;
	return $self->{_blast};
}

sub setBLASTstatus {
	my ($self, $blast) = @_;
 	$self->{_blast} = $blast;
}

sub setRelationToPhage {
	my ($self, $relation) = @_;
 	$self->{_namerelated} = $relation;
}

sub isRelatedToPhage {
	my ($self) = @_;
	return $self->{_namerelated};
}

sub getStart {
	my ($self) = @_;
	return $self->{_start};
}

sub getEnd {
	my ($self) = @_;
	return $self->{_end};
}

sub getLocalGI {
	my ($self) = @_;
	return $self->{_localgi};
}

sub getProduct {
	my ($self) = @_;
	return $self->{_product};
}

sub isRNA {
	my ($self) = @_;
	if ($self->{_type} eq 't'){
		return 1;
	}
	elsif ($self->{_type} eq 'tm'){
		return 2;
	}
	elsif ($self->{_type} eq 'r'){
		return 3;
	}
	return 0;
}

sub is5to3 {
	my ($self) = @_;
	if ($self->{_strand} eq '+'){
		return 1;
	}
	return 0;
}

sub next {
	my ($self) = @_;
	return $self->{_next};
}

sub previous {
	my ($self) = @_;
	return $self->{_previous};
}

sub setNext {
	my ( $self, $next ) = @_;
    $self->{_next} = $next;
}

sub setPrevious {
	my ( $self, $previous ) = @_;
    $self->{_previous} = $previous;
}

sub printInfo2 {
	my ( $self, $next ) = @_;
	print "yes $next end\n";
}

sub isInGene {
	my ($self, $pos) = @_;
	if ($pos >= $self->getStart() and $pos <= $self->getEnd()){
		return 1;
	}
	return 0;
}

# Get contig id, indicating contig to which this gene belongs.
sub getContigId {
	my ($self) = @_;
	return $self->{_contigId};
}

# Set contig id, indicating contig to which this gene belongs.
sub setContigId {
	my ( $self, $contigId ) = @_;
    $self->{_contigId} = $contigId;
}

sub printInfo {
	my ($self, $flag) =  @_;
	if ($self->{_hit} eq ''){
		print ' ';
	}
	else{
		print '*'
	}
	if ($self->{_namerelated} == 0){
		print ' ';
	}
	else{
		print 'x'
	}
	if ($self->isIntegrase()){
		print 'i';
	}
	elsif ($self->isTransposase()){
		print 't';
	}
	else{
		print ' ';
	}
	if ($self->isRNA() == 1){
		print "  tRNA  ";
	}
	elsif ($self->isRNA() == 2){
		print "  tmRNA ";
	}
	elsif ($self->isRNA() == 3){
		print "  rRNA  ";
	}
	else{
		print "  prot  "
	}
	if ($self->is5to3()){
		print $self->getStart(), ",", $self->getEnd(), " ";
	}
	else{
		print "complement(", $self->getStart(), ",", $self->getEnd(),") ";
	}
	print $self->{_localgi},"  ",$self->getProduct(),"\n";
	if ($flag == 1){
		my $hits = $self->getBLASTresult();
		if ($hits ne '' and $hits ne "keyword"){
			foreach (@$hits){
				$_->printInfo();
			}
		}
	}
}

sub isIntegrase {
	my ($self) = @_;
	if ($self->getProduct() =~m/integrase|specific recombinase|phage recombinase/i){
		return 1;
	}
	my $bla = $self->getBLASTresult();
	if ($bla ne '' and $bla ne 'keyword'){
		foreach (@$bla){
			if ($_->isIntegrase() == 1){
				return 1;
			}
		}
	}
	return 0;
}

sub isTransposase {
	my ($self) = @_;
	my $bla = $self->getBLASTresult();
	if ($bla eq "keyword"){
		if ($self->getProduct() =~m/transposase/i){
			return 1;
		}
	}
	elsif ($bla ne ''){
		foreach (@$bla){
			if ($_->isTransposase() == 1){
				return 1;
			}
		}
	}
	return 0;
}

1;
