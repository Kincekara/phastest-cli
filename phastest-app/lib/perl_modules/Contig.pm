#!/usr/bin/perl

# David Arndt, April 2016
# A class to represent a contig.

use strict;
package Contig;

# build a contig
sub new
{
	my $class = shift;
	my $self = {
		_start => shift, 	# Start position within concatenated contigs sequence (inclusive).
			# First nucleotide in concatenated sequence has position 1.
		_end => shift, 		# End position within concatenated contigs sequence (inclusive).
		_label => shift,	# Label given to the contig, including starting > character.
		_contigId => -1 # id used to identify the contig
	};
	bless $self, $class;
	return $self;
}

sub getStart {
	my ($self) = @_;
	return $self->{_start};
}

sub getEnd {
	my ($self) = @_;
	return $self->{_end};
}

sub getLabel {
	my ($self) = @_;
	return $self->{_label};
}

sub getId {
	my ($self) = @_;
	return $self->{_contigId};
}

sub setId {
	my ( $self, $id ) = @_;
    $self->{_contigId} = $id;
}

sub isInContig {
	my ($self, $pos) = @_;
	if ($pos >= $self->getStart() and $pos <= $self->getEnd()){
		return 1;
	}
	return 0;
}

1;
