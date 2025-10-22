#!/usr/bin/perl
package AttachmentSite;
sub new
{
	my $class = shift;
    my $self = {
		_atts => '',	# array of attachment sites
		_attL => '',
		_attR => '',
		_length => 0
	};
	bless $self, $class;
    return $self;
}

 
1;