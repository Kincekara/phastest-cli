#!/usr/bin/perl
# April 13, 2011 by Rah: handle the sequence file
package Status;
use strict;

sub new
{
	my ($class) = shift;
	my @array = ();
	my $self = {
		_prophage => shift,
		_int_tran => shift,
		_score => 0,
		_att => \@array
	};
	bless $self, $class;
    return $self;
}

sub cp {
	my ($self, $to) = @_;
	$to->{_score} = $self->{_score};
	foreach (@{$self->{_att}}){
		push @{$to->{_att}}, $_;
	}
}
# atts must fully include or exclude each other and no more att
# than int+transposase 
sub add {
	my ($self, $att) = @_;
#	print "add $att->{_attL} $att->{_attR} ";
	my $data = $self->{_att};
	my $intandtran = $self->{_int_tran};
	if ($#$data >= $#$intandtran){
		return -1;
	}
	foreach (@$data){
		if ($att->{_attL} < $_->{_attL} and $att->{_attR} > $_->{_attL} and $att->{_attR} < $_->{_attR}){
#			print " BAD\n";
			return -1;
		}
		if ($att->{_attL} > $_->{_attL} and $att->{_attL} < $_->{_attR} and $att->{_attR} > $_->{_attR}){
#			print " BAD\n";
			return -1;
		}
	}
#	print " OK\n";
	push @$data, $att;
	if ($self->hasContents() == 0){
		pop (@$data);
		return -1;
	}
	return 1;
}

sub rm {
	my ($self, $att) = @_;
	my $array = $self->{_att};
	if ($$array[$#$array] == $att){
		pop (@$array);
	}
	else{
		die "try to rm unexpected value from stack";
	}
}

sub getCombination {
	my ($self) = @_;
	return $self->{_att};
}

# make sure every att has include at least one int 
sub hasContents{
	my ($self) = @_;
	my @atts = ();
	my @cores = ();
	push @atts, @{$self->{_att}};
	push @cores, @{$self->{_int_tran}};
	# pick up innerest att pair
	while (&rmInnerestAtt(\@atts, \@cores) == 1){};
#	if ($#atts>-1){exit;}
	foreach (@atts){
		if ($_ ne ''){
			return 0;
		}
	}
	return 1;
}

sub rmInnerestAtt {
	my ($att, $cores) = @_;
#	print "att has $#$att+1   core has $#$cores+1\n";
	for my $i (0 .. $#$att){
		if ($$att[$i] eq ''){
			next;
		}
		my $innerest = 1;
		foreach (@$att){
			if ($_ eq '' or $_ == $$att[$i]){
				next;
			}
			if ($_->{_attL} >= $$att[$i]->{_attL} and $_->{_attR} <= $$att[$i]->{_attR}){
				$innerest = 0;
				last;
			}
		}
		if ($innerest == 0){
			next;
		}
		my $hasCore = 0;
		for my $n (0 .. $#$cores){
			if ($$cores[$n] eq ''){
				next;
			}
			if ($$att[$i]->{_attL} < $$cores[$n]->getStart() and $$att[$i]->{_attR} > $$cores[$n]->getEnd()){
				$hasCore = 1;
				$$cores[$n] = '';
			}
		}
		if ($hasCore == 1){
			$$att[$i] = '';
			return 1;
		}
		return 0;
	}
	return 0;
}
sub printInfo {
	my ($self) = @_;
	my $array = $self->{_att};
	foreach (@$array){
		print "($_->{_attL} , $_->{_attR}) ";
	}
	print "<",$self->hasContents(),",",$self->getScore(),">\n";
}

# get score based on how many phage gene included
sub getScore {
	my ($self) = @_;
	my $score = 0;
	for (my $cursor = $self->{_prophage}->getHead(); $cursor ne ''; $cursor = $cursor->next()){
		foreach (@{$self->{_att}}){
			if ($_->{_attL} < $cursor->getStart() and $_->{_attR} > $cursor->getEnd()){
				$score++;
			}
		}
		if ($cursor == $self->{_prophage}->getTail()){
			last;
		}
	}
	return $score;
}
1;
