package BiochemistryAPI::BiochemistryAPIClient;

use JSON::RPC::Client;
use POSIX;
use strict;
use Data::Dumper;
use URI;
use Bio::KBase::Exceptions;
my $get_time = sub { time, 0 };
eval {
    require Time::HiRes;
    $get_time = sub { Time::HiRes::gettimeofday() };
};

use Bio::KBase::AuthToken;

# Client version should match Impl version
# This is a Semantic Version number,
# http://semver.org
our $VERSION = "0.1.0";

=head1 NAME

BiochemistryAPI::BiochemistryAPIClient

=head1 DESCRIPTION


A KBase module: BiochemistryAPI


=cut

sub new
{
    my($class, $url, @args) = @_;
    

    my $self = {
	client => BiochemistryAPI::BiochemistryAPIClient::RpcClient->new,
	url => $url,
	headers => [],
    };

    chomp($self->{hostname} = `hostname`);
    $self->{hostname} ||= 'unknown-host';

    #
    # Set up for propagating KBRPC_TAG and KBRPC_METADATA environment variables through
    # to invoked services. If these values are not set, we create a new tag
    # and a metadata field with basic information about the invoking script.
    #
    if ($ENV{KBRPC_TAG})
    {
	$self->{kbrpc_tag} = $ENV{KBRPC_TAG};
    }
    else
    {
	my ($t, $us) = &$get_time();
	$us = sprintf("%06d", $us);
	my $ts = strftime("%Y-%m-%dT%H:%M:%S.${us}Z", gmtime $t);
	$self->{kbrpc_tag} = "C:$0:$self->{hostname}:$$:$ts";
    }
    push(@{$self->{headers}}, 'Kbrpc-Tag', $self->{kbrpc_tag});

    if ($ENV{KBRPC_METADATA})
    {
	$self->{kbrpc_metadata} = $ENV{KBRPC_METADATA};
	push(@{$self->{headers}}, 'Kbrpc-Metadata', $self->{kbrpc_metadata});
    }

    if ($ENV{KBRPC_ERROR_DEST})
    {
	$self->{kbrpc_error_dest} = $ENV{KBRPC_ERROR_DEST};
	push(@{$self->{headers}}, 'Kbrpc-Errordest', $self->{kbrpc_error_dest});
    }

    #
    # This module requires authentication.
    #
    # We create an auth token, passing through the arguments that we were (hopefully) given.

    {
	my %arg_hash2 = @args;
	if (exists $arg_hash2{"token"}) {
	    $self->{token} = $arg_hash2{"token"};
	} elsif (exists $arg_hash2{"user_id"}) {
	    my $token = Bio::KBase::AuthToken->new(@args);
	    if (!$token->error_message) {
	        $self->{token} = $token->token;
	    }
	}
	
	if (exists $self->{token})
	{
	    $self->{client}->{token} = $self->{token};
	}
    }

    my $ua = $self->{client}->ua;	 
    my $timeout = $ENV{CDMI_TIMEOUT} || (30 * 60);	 
    $ua->timeout($timeout);
    bless $self, $class;
    #    $self->_validate_version();
    return $self;
}




=head2 get_reactions

  $out_reactions = $obj->get_reactions($params)

=over 4

=item Parameter and return types

=begin html

<pre>
$params is a BiochemistryAPI.get_reactions_params
$out_reactions is a reference to a list where each element is a BiochemistryAPI.Reaction
get_reactions_params is a reference to a hash where the following keys are defined:
	reactions has a value which is a reference to a list where each element is a BiochemistryAPI.reaction_id
reaction_id is a string
Reaction is a reference to a hash where the following keys are defined:
	id has a value which is a BiochemistryAPI.reaction_id
	name has a value which is a string
	abbrev has a value which is a string
	enzymes has a value which is a reference to a list where each element is a string
	direction has a value which is a string
	reversibility has a value which is a string
	deltaG has a value which is a float
	deltaGErr has a value which is a float
	equation has a value which is a string
	definition has a value which is a string

</pre>

=end html

=begin text

$params is a BiochemistryAPI.get_reactions_params
$out_reactions is a reference to a list where each element is a BiochemistryAPI.Reaction
get_reactions_params is a reference to a hash where the following keys are defined:
	reactions has a value which is a reference to a list where each element is a BiochemistryAPI.reaction_id
reaction_id is a string
Reaction is a reference to a hash where the following keys are defined:
	id has a value which is a BiochemistryAPI.reaction_id
	name has a value which is a string
	abbrev has a value which is a string
	enzymes has a value which is a reference to a list where each element is a string
	direction has a value which is a string
	reversibility has a value which is a string
	deltaG has a value which is a float
	deltaGErr has a value which is a float
	equation has a value which is a string
	definition has a value which is a string


=end text

=item Description

Returns data for the requested reactions

=back

=cut

 sub get_reactions
{
    my($self, @args) = @_;

# Authentication: required

    if ((my $n = @args) != 1)
    {
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
							       "Invalid argument count for function get_reactions (received $n, expecting 1)");
    }
    {
	my($params) = @args;

	my @_bad_arguments;
        (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"params\" (value was \"$params\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to get_reactions:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
								   method_name => 'get_reactions');
	}
    }

    my $url = $self->{url};
    my $result = $self->{client}->call($url, $self->{headers}, {
	    method => "BiochemistryAPI.get_reactions",
	    params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
					       code => $result->content->{error}->{code},
					       method_name => 'get_reactions',
					       data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
					      );
	} else {
	    return wantarray ? @{$result->result} : $result->result->[0];
	}
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method get_reactions",
					    status_line => $self->{client}->status_line,
					    method_name => 'get_reactions',
				       );
    }
}
 


=head2 get_compounds

  $out_compounds = $obj->get_compounds($params)

=over 4

=item Parameter and return types

=begin html

<pre>
$params is a BiochemistryAPI.get_compounds_params
$out_compounds is a reference to a list where each element is a BiochemistryAPI.Compound
get_compounds_params is a reference to a hash where the following keys are defined:
	compounds has a value which is a reference to a list where each element is a BiochemistryAPI.compound_id
compound_id is a string
Compound is a reference to a hash where the following keys are defined:
	id has a value which is a BiochemistryAPI.compound_id
	abbrev has a value which is a string
	name has a value which is a string
	aliases has a value which is a reference to a list where each element is a string
	charge has a value which is a float
	deltaG has a value which is a float
	deltaGErr has a value which is a float
	formula has a value which is a string

</pre>

=end html

=begin text

$params is a BiochemistryAPI.get_compounds_params
$out_compounds is a reference to a list where each element is a BiochemistryAPI.Compound
get_compounds_params is a reference to a hash where the following keys are defined:
	compounds has a value which is a reference to a list where each element is a BiochemistryAPI.compound_id
compound_id is a string
Compound is a reference to a hash where the following keys are defined:
	id has a value which is a BiochemistryAPI.compound_id
	abbrev has a value which is a string
	name has a value which is a string
	aliases has a value which is a reference to a list where each element is a string
	charge has a value which is a float
	deltaG has a value which is a float
	deltaGErr has a value which is a float
	formula has a value which is a string


=end text

=item Description

Returns data for the requested compounds

=back

=cut

 sub get_compounds
{
    my($self, @args) = @_;

# Authentication: required

    if ((my $n = @args) != 1)
    {
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
							       "Invalid argument count for function get_compounds (received $n, expecting 1)");
    }
    {
	my($params) = @args;

	my @_bad_arguments;
        (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"params\" (value was \"$params\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to get_compounds:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
								   method_name => 'get_compounds');
	}
    }

    my $url = $self->{url};
    my $result = $self->{client}->call($url, $self->{headers}, {
	    method => "BiochemistryAPI.get_compounds",
	    params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
					       code => $result->content->{error}->{code},
					       method_name => 'get_compounds',
					       data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
					      );
	} else {
	    return wantarray ? @{$result->result} : $result->result->[0];
	}
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method get_compounds",
					    status_line => $self->{client}->status_line,
					    method_name => 'get_compounds',
				       );
    }
}
 


=head2 substructure_search

  $matching_ids = $obj->substructure_search($params)

=over 4

=item Parameter and return types

=begin html

<pre>
$params is a BiochemistryAPI.substructure_search_params
$matching_ids is a reference to a list where each element is a BiochemistryAPI.compound_id
substructure_search_params is a reference to a hash where the following keys are defined:
	query has a value which is a string
compound_id is a string

</pre>

=end html

=begin text

$params is a BiochemistryAPI.substructure_search_params
$matching_ids is a reference to a list where each element is a BiochemistryAPI.compound_id
substructure_search_params is a reference to a hash where the following keys are defined:
	query has a value which is a string
compound_id is a string


=end text

=item Description

Returns compound ids for compounds that contain the query substructure

=back

=cut

 sub substructure_search
{
    my($self, @args) = @_;

# Authentication: required

    if ((my $n = @args) != 1)
    {
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
							       "Invalid argument count for function substructure_search (received $n, expecting 1)");
    }
    {
	my($params) = @args;

	my @_bad_arguments;
        (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"params\" (value was \"$params\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to substructure_search:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
								   method_name => 'substructure_search');
	}
    }

    my $url = $self->{url};
    my $result = $self->{client}->call($url, $self->{headers}, {
	    method => "BiochemistryAPI.substructure_search",
	    params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
					       code => $result->content->{error}->{code},
					       method_name => 'substructure_search',
					       data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
					      );
	} else {
	    return wantarray ? @{$result->result} : $result->result->[0];
	}
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method substructure_search",
					    status_line => $self->{client}->status_line,
					    method_name => 'substructure_search',
				       );
    }
}
 


=head2 similarity_search

  $matching_ids = $obj->similarity_search($params)

=over 4

=item Parameter and return types

=begin html

<pre>
$params is a BiochemistryAPI.similarity_search_params
$matching_ids is a reference to a list where each element is a BiochemistryAPI.compound_id
similarity_search_params is a reference to a hash where the following keys are defined:
	query has a value which is a string
	fp_type has a value which is a string
	min_similarity has a value which is a float
compound_id is a string

</pre>

=end html

=begin text

$params is a BiochemistryAPI.similarity_search_params
$matching_ids is a reference to a list where each element is a BiochemistryAPI.compound_id
similarity_search_params is a reference to a hash where the following keys are defined:
	query has a value which is a string
	fp_type has a value which is a string
	min_similarity has a value which is a float
compound_id is a string


=end text

=item Description

Returns compound ids for compounds that have greater fingerprint similarity than the min_similarity threshold

=back

=cut

 sub similarity_search
{
    my($self, @args) = @_;

# Authentication: required

    if ((my $n = @args) != 1)
    {
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
							       "Invalid argument count for function similarity_search (received $n, expecting 1)");
    }
    {
	my($params) = @args;

	my @_bad_arguments;
        (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"params\" (value was \"$params\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to similarity_search:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
								   method_name => 'similarity_search');
	}
    }

    my $url = $self->{url};
    my $result = $self->{client}->call($url, $self->{headers}, {
	    method => "BiochemistryAPI.similarity_search",
	    params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
					       code => $result->content->{error}->{code},
					       method_name => 'similarity_search',
					       data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
					      );
	} else {
	    return wantarray ? @{$result->result} : $result->result->[0];
	}
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method similarity_search",
					    status_line => $self->{client}->status_line,
					    method_name => 'similarity_search',
				       );
    }
}
 


=head2 depict_compounds

  $depictions = $obj->depict_compounds($params)

=over 4

=item Parameter and return types

=begin html

<pre>
$params is a BiochemistryAPI.depict_compounds_params
$depictions is a reference to a list where each element is a string
depict_compounds_params is a reference to a hash where the following keys are defined:
	compound_structures has a value which is a reference to a list where each element is a string

</pre>

=end html

=begin text

$params is a BiochemistryAPI.depict_compounds_params
$depictions is a reference to a list where each element is a string
depict_compounds_params is a reference to a hash where the following keys are defined:
	compound_structures has a value which is a reference to a list where each element is a string


=end text

=item Description

Returns a list of depictions for the compound_structures in SVG format

=back

=cut

 sub depict_compounds
{
    my($self, @args) = @_;

# Authentication: required

    if ((my $n = @args) != 1)
    {
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
							       "Invalid argument count for function depict_compounds (received $n, expecting 1)");
    }
    {
	my($params) = @args;

	my @_bad_arguments;
        (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"params\" (value was \"$params\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to depict_compounds:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
								   method_name => 'depict_compounds');
	}
    }

    my $url = $self->{url};
    my $result = $self->{client}->call($url, $self->{headers}, {
	    method => "BiochemistryAPI.depict_compounds",
	    params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
					       code => $result->content->{error}->{code},
					       method_name => 'depict_compounds',
					       data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
					      );
	} else {
	    return wantarray ? @{$result->result} : $result->result->[0];
	}
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method depict_compounds",
					    status_line => $self->{client}->status_line,
					    method_name => 'depict_compounds',
				       );
    }
}
 
  
sub status
{
    my($self, @args) = @_;
    if ((my $n = @args) != 0) {
        Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
                                   "Invalid argument count for function status (received $n, expecting 0)");
    }
    my $url = $self->{url};
    my $result = $self->{client}->call($url, $self->{headers}, {
        method => "BiochemistryAPI.status",
        params => \@args,
    });
    if ($result) {
        if ($result->is_error) {
            Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
                           code => $result->content->{error}->{code},
                           method_name => 'status',
                           data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
                          );
        } else {
            return wantarray ? @{$result->result} : $result->result->[0];
        }
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method status",
                        status_line => $self->{client}->status_line,
                        method_name => 'status',
                       );
    }
}
   

sub version {
    my ($self) = @_;
    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
        method => "BiochemistryAPI.version",
        params => [],
    });
    if ($result) {
        if ($result->is_error) {
            Bio::KBase::Exceptions::JSONRPC->throw(
                error => $result->error_message,
                code => $result->content->{code},
                method_name => 'depict_compounds',
            );
        } else {
            return wantarray ? @{$result->result} : $result->result->[0];
        }
    } else {
        Bio::KBase::Exceptions::HTTP->throw(
            error => "Error invoking method depict_compounds",
            status_line => $self->{client}->status_line,
            method_name => 'depict_compounds',
        );
    }
}

sub _validate_version {
    my ($self) = @_;
    my $svr_version = $self->version();
    my $client_version = $VERSION;
    my ($cMajor, $cMinor) = split(/\./, $client_version);
    my ($sMajor, $sMinor) = split(/\./, $svr_version);
    if ($sMajor != $cMajor) {
        Bio::KBase::Exceptions::ClientServerIncompatible->throw(
            error => "Major version numbers differ.",
            server_version => $svr_version,
            client_version => $client_version
        );
    }
    if ($sMinor < $cMinor) {
        Bio::KBase::Exceptions::ClientServerIncompatible->throw(
            error => "Client minor version greater than Server minor version.",
            server_version => $svr_version,
            client_version => $client_version
        );
    }
    if ($sMinor > $cMinor) {
        warn "New client version available for BiochemistryAPI::BiochemistryAPIClient\n";
    }
    if ($sMajor == 0) {
        warn "BiochemistryAPI::BiochemistryAPIClient version is $svr_version. API subject to change.\n";
    }
}

=head1 TYPES



=head2 compound_id

=over 4



=item Description

An identifier for compounds in the KBase biochemistry database. e.g. cpd00001


=item Definition

=begin html

<pre>
a string
</pre>

=end html

=begin text

a string

=end text

=back



=head2 reaction_id

=over 4



=item Description

A string identifier used for a reaction in a KBase biochemistry.


=item Definition

=begin html

<pre>
a string
</pre>

=end html

=begin text

a string

=end text

=back



=head2 Compound

=over 4



=item Description

Data structures for media formulation

                compound_id id - ID of compound
                string abbrev - abbreviated name of compound
                string name - primary name of compound
                list<string> aliases - list of aliases for compound
                float charge - molecular charge of compound
                float deltaG - estimated compound delta G
                float deltaGErr - uncertainty in estimated compound delta G
                string formula - molecular formula of compound


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
id has a value which is a BiochemistryAPI.compound_id
abbrev has a value which is a string
name has a value which is a string
aliases has a value which is a reference to a list where each element is a string
charge has a value which is a float
deltaG has a value which is a float
deltaGErr has a value which is a float
formula has a value which is a string

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
id has a value which is a BiochemistryAPI.compound_id
abbrev has a value which is a string
name has a value which is a string
aliases has a value which is a reference to a list where each element is a string
charge has a value which is a float
deltaG has a value which is a float
deltaGErr has a value which is a float
formula has a value which is a string


=end text

=back



=head2 Reaction

=over 4



=item Description

Data structures for media formulation

                reaction_id id - ID of reaction
                string name - primary name of reaction
                string abbrev - abbreviated name of reaction
                list<string> enzymes - list of EC numbers for reaction
                string direction - directionality of reaction
                string reversibility - reversibility of reaction
                float deltaG - estimated delta G of reaction
                float deltaGErr - uncertainty in estimated delta G of reaction
                string equation - reaction equation in terms of compound IDs
                string definition - reaction equation in terms of compound names


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
id has a value which is a BiochemistryAPI.reaction_id
name has a value which is a string
abbrev has a value which is a string
enzymes has a value which is a reference to a list where each element is a string
direction has a value which is a string
reversibility has a value which is a string
deltaG has a value which is a float
deltaGErr has a value which is a float
equation has a value which is a string
definition has a value which is a string

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
id has a value which is a BiochemistryAPI.reaction_id
name has a value which is a string
abbrev has a value which is a string
enzymes has a value which is a reference to a list where each element is a string
direction has a value which is a string
reversibility has a value which is a string
deltaG has a value which is a float
deltaGErr has a value which is a float
equation has a value which is a string
definition has a value which is a string


=end text

=back



=head2 get_reactions_params

=over 4



=item Description

Input parameters for the "get_reactions" function.

                list<reaction_id> reactions - a list of the reaction IDs for the reactions to be returned (a required argument)


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
reactions has a value which is a reference to a list where each element is a BiochemistryAPI.reaction_id

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
reactions has a value which is a reference to a list where each element is a BiochemistryAPI.reaction_id


=end text

=back



=head2 get_compounds_params

=over 4



=item Description

Input parameters for the "get_compounds" function.
    list<compound_id> compounds - a list of the compound IDs for the compounds to be returned (a required argument)


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
compounds has a value which is a reference to a list where each element is a BiochemistryAPI.compound_id

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
compounds has a value which is a reference to a list where each element is a BiochemistryAPI.compound_id


=end text

=back



=head2 substructure_search_params

=over 4



=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
query has a value which is a string

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
query has a value which is a string


=end text

=back



=head2 similarity_search_params

=over 4



=item Description

string query: Either InChI or SMILES string
        string fp_type: Either MACCS or Morgan fingerprints
        float min_similarity: In range 0-1


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
query has a value which is a string
fp_type has a value which is a string
min_similarity has a value which is a float

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
query has a value which is a string
fp_type has a value which is a string
min_similarity has a value which is a float


=end text

=back



=head2 depict_compounds_params

=over 4



=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
compound_structures has a value which is a reference to a list where each element is a string

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
compound_structures has a value which is a reference to a list where each element is a string


=end text

=back



=cut

package BiochemistryAPI::BiochemistryAPIClient::RpcClient;
use base 'JSON::RPC::Client';
use POSIX;
use strict;

#
# Override JSON::RPC::Client::call because it doesn't handle error returns properly.
#

sub call {
    my ($self, $uri, $headers, $obj) = @_;
    my $result;


    {
	if ($uri =~ /\?/) {
	    $result = $self->_get($uri);
	}
	else {
	    Carp::croak "not hashref." unless (ref $obj eq 'HASH');
	    $result = $self->_post($uri, $headers, $obj);
	}

    }

    my $service = $obj->{method} =~ /^system\./ if ( $obj );

    $self->status_line($result->status_line);

    if ($result->is_success) {

        return unless($result->content); # notification?

        if ($service) {
            return JSON::RPC::ServiceObject->new($result, $self->json);
        }

        return JSON::RPC::ReturnObject->new($result, $self->json);
    }
    elsif ($result->content_type eq 'application/json')
    {
        return JSON::RPC::ReturnObject->new($result, $self->json);
    }
    else {
        return;
    }
}


sub _post {
    my ($self, $uri, $headers, $obj) = @_;
    my $json = $self->json;

    $obj->{version} ||= $self->{version} || '1.1';

    if ($obj->{version} eq '1.0') {
        delete $obj->{version};
        if (exists $obj->{id}) {
            $self->id($obj->{id}) if ($obj->{id}); # if undef, it is notification.
        }
        else {
            $obj->{id} = $self->id || ($self->id('JSON::RPC::Client'));
        }
    }
    else {
        # $obj->{id} = $self->id if (defined $self->id);
	# Assign a random number to the id if one hasn't been set
	$obj->{id} = (defined $self->id) ? $self->id : substr(rand(),2);
    }

    my $content = $json->encode($obj);

    $self->ua->post(
        $uri,
        Content_Type   => $self->{content_type},
        Content        => $content,
        Accept         => 'application/json',
	@$headers,
	($self->{token} ? (Authorization => $self->{token}) : ()),
    );
}



1;
