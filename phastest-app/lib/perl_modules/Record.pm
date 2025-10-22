#!/usr/bin/perl
# April 12, 2011 by Rah: a class of blast returned annotation
use strict;
package Record;
sub new
{
	my $class = shift;
	my $self = {
		# permental data
		_gi=>shift,
		_species => shift, 	# phage species
		_definition => shift, 	# definition in database
		_ref => shift, # ref acc #
		_evalue => shift,	# evalue
		_querygi => shift       # gi number of hit
	};
	bless $self, $class;
	return $self;
}
sub getQueryGI{
	my $self=shift;
	return $self->{_querygi};
}

sub getEvalue {
	my ($self) = @_;
	return $self->{_evalue};
}

sub getRefAcc {
	my ($self) = @_;
	return $self->{_ref};
}
sub setSpecies {
	my ($self, $spe) = @_;
	$self->{_species} = $spe;
}

sub setDefinition{
	my ($self, $def) = @_;
	$self->{_definition} = $def;
}

sub getSpecies {
	my ($self) = @_;
	return $self->{_species};
}

sub getDefinition{
	my ($self) = @_;
	return $self->{_definition};
}

sub getGI {
	my ($self) = @_;
	return $self->{_gi};
}
sub printInfo {
	my ($self) = @_;
	print "        ACC:", $self->getRefAcc(), " SPE:", $self->getSpecies(), " DEF:", $self->getDefinition(), "\n";
}

sub isIntegrase {
	my ($self) = @_;
	if ($self->{_definition} =~m/integrase/i or $self->{_definition} =~m/site.specific recombinase/i or $self->{_definition} =~m/phage recombinase/i){
		return 1;
	}
	return 0;
}

sub isTransposase {
	my ($self) = @_;
	if ($self->{_definition} =~m/transposase/i){
		return 1;
	}
	return 0;
}
1;
