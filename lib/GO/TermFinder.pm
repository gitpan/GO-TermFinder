package GO::TermFinder;

# File        : TermFinder.pm
# Author      : Gavin Sherlock
# Date Begun  : December 31st 2002

# $Id: TermFinder.pm,v 1.16 2003/03/03 16:50:24 sherlock Exp $

# License information (the MIT license)

# Copyright (c) 2003 Gavin Sherlock; Stanford University

# Permission is hereby granted, free of charge, to any person
# obtaining a copy of this software and associated documentation files
# (the "Software"), to deal in the Software without restriction,
# including without limitation the rights to use, copy, modify, merge,
# publish, distribute, sublicense, and/or sell copies of the Software,
# and to permit persons to whom the Software is furnished to do so,
# subject to the following conditions:

# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT.  IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
# BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
# ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
# CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

=pod

=head1 NAME

GO::TermFinder

=head1 DESCRIPTION

This package is intended to provide a method whereby the P-values of a
set of GO annotations can be determined for a set of genes, based on
the number of genes that exist in the particular genome, and their
annotation, and the frequency with which the GO nodes are annotated
across the provided set of genes.  The P-value is simply calculated
using either the hypergeometric, or the binomial distribution, as the
probability of x or more out of n genes having a given annotation,
given that G of N have that annotation in the genome in general.  The
hypergeometric distribution (sampling without replacement) is more
accurate, though slower to calculate the the binomial distibution
(sampling with replacement).

In addition, a corrected p-value is also calculated, to correct for
multiple hypothesis testing.  From the sub-graph of GO nodes that are
considered as hypotheses (any nodes with 2 or more annotations from
the list of genes to be considered), the minimal set of nodes is
determined from which all other hypotheses (nodes and annotations) can
be inferred.  Their number is then used to multiply the calculated
p-values, to generate corrected p-values.  The client has access to
both the corrected an uncorrected values.  This is done instead of
simple Bonferroni correction, because the hypotheses are not
independent of one another.

The general idea is that a list of genes may have been identified for
some reason, e.g. they are coregulated, and TermFinder can be used to
found out if any nodes annotate the set of genes to a level which is
extremely improbable if the genes had simply been picked at random.

=head1 TODO

On the day that I finished the code, BioPerl 1.2 was released, which
appears to have some Ontology parsing code that is better than mine,
which uses the Graph module on CPAN to implement the DAG structure of
GO.  The Graph module provides a richer set of methods for querying
the structure, though the BioPerl implementation does not appear to
have all the things I want.  Anyway, converting my code to use the
BioPerl Ontology stuff is probably a good goal, as the modules in
their appear to have been somewhat better conceived in terms of their
relationships to each other.  I did not notice anything to do with
providing annotations, so my
GO::AnnotationProvider::AnnotationParser may still be of use.  We
will see....

May want the client to decide the behaviour for ambiguous names,
rather than having it hard coded (eg always ignore; use if standard
name (current implementation); use all databaseIds for the ambiguous
name; decide on a case by case basis (potentially useful if running on
command line)).

=cut

use strict;
use warnings;
use diagnostics;

use vars qw ($PACKAGE $VERSION);

use GO::Node;

$VERSION = 0.1;
$PACKAGE = 'GO:TermFinder';

# class variables

my @kRequiredArgs = qw (annotationProvider ontologyProvider totalNumGenes aspect);

my $kArgs                   = $PACKAGE.'::__args';
my $kTotalGoNodeCounts      = $PACKAGE.'::__totalGoNodeCounts';
my $kGoCounts               = $PACKAGE.'::__goCounts';
my $kGOIDsForDatabaseIds    = $PACKAGE.'::__goidsForDatabaseIds';
my $kDatabaseIds            = $PACKAGE.'::__databaseIds';
my $kTotalNumAnnotatedGenes = $PACKAGE.'::__totalNumAnnotatedGenes';
my $kMethod                 = $PACKAGE.'::__method';
my $kLogFactorials          = $PACKAGE.'::__logFactorials';
my $kLogNCr                 = $PACKAGE.'::__logNCr';
my $kPvalues                = $PACKAGE.'::__pValues';

my %kAllowedMethods = ('hypergeometric' => undef,
		       'binomial'       => undef); # the methods by which the p-value can be calculated

# set up a GO node that corresponds to anything passed in that has no
# annotation

my $kUnannotatedNode = GO::Node->new(goid => "GO:XXXXXXX",
				     term => "unannotated");

#####################################################################
sub new{
#####################################################################
# This is the constructor.  It expects to be passed named arguments
# for an annotationProvider, an ontologyProvider, and how many genes
# in total exist.  In addition, it must be told the aspect of the ontology
# provider, so that it knows how to query the annotationProvider.
#
# Usage :
#
#    my $termFinder = GO::TermFinder->new(annotationProvider=> $annotationProvider,
#                                         ontologyProvider  => $ontologyProvider,
#                                         totalNumGenes     => $num,
#                                         aspect            => <P|C|F>);
#

    my ($class, %args) = @_;

    my $self = {};

    bless $self, $class;

    $self->__checkAndStoreArgs(%args);

    $self->__init; # initialize counts for all GO nodes

    return $self;

}

#####################################################################
sub __checkAndStoreArgs{
#####################################################################
# This private method simply checks that all the required arguments
# have been provided, and stores them within the object

    my ($self, %args) = @_;

    foreach my $arg (@kRequiredArgs){

	if (!exists ($args{$arg})){

	    die "You did not provide a $arg argument.";

	}elsif (!defined ($args{$arg})){

	    die "Your $arg argument is not defined";

	}

	$self->{$kArgs}{$arg} = $args{$arg}; # store in object

    }   

}

#####################################################################
sub __init{
#####################################################################
# This private method determines all counts to all GO nodes, as the
# background frequency of annotations in the genome

    my ($self) = @_;

    my @allDatabaseIds = $self->__annotationProvider->allDatabaseIds();

    my $totalNumAnnotatedGenes = scalar(@allDatabaseIds);

    # check that they said there's at least as many genes in total
    # as the annotation provider says that there is.    

    if ($totalNumAnnotatedGenes > $self->__totalNumGenes){

	print "The annotation provider indicates that there are more genes than the client indicated.\n";
	print "The annotaion provider indicates there are $totalNumAnnotatedGenes, while the client indicated only ", $self->__totalNumGenes, ".\n";
	print "Thus assuming the total number of genes is that indicated by the annotation provider.\n";

	$self->{$kArgs}{totalNumGenes} = $totalNumAnnotatedGenes;

    }

    my $totalNodeCounts = $self->__buildHashRefOfAnnotations(\@allDatabaseIds);

    if ($totalNumAnnotatedGenes < $self->__totalNumGenes){

    	# if there are extra, entirely unannotated genes (indicated by
    	# the total number of genes provided being greater than the
    	# number that existed in the annotation provider), we must
    	# make sure that it's treated that they will at least be
    	# annotated to the root (Gene Ontology), and its immediate
    	# child (which is the name of the Ontology, eg
    	# Biological_process, Molecular_function, and
    	# Cellular_component), and the 'unannotated' node

	# so simply add extra annotations

	my $rootNodeId  = $self->__ontologyProvider->rootNode->goid;
	
	my $childNodeId = ($self->__ontologyProvider->rootNode->childNodes())[0]->goid;

	$totalNodeCounts->{$rootNodeId} = $self->__totalNumGenes;

	$totalNodeCounts->{$childNodeId} += ($self->__totalNumGenes - $totalNumAnnotatedGenes);

	$totalNodeCounts->{$kUnannotatedNode->goid} += ($self->__totalNumGenes - $totalNumAnnotatedGenes);

    }

    $self->{$kTotalGoNodeCounts}      = $totalNodeCounts;
    $self->{$kTotalNumAnnotatedGenes} = $totalNumAnnotatedGenes;

    $self->__cacheLogFactorials;

}

#####################################################################
sub __cacheLogFactorials{
#####################################################################
# The maximum factorial that will ever have to be calculated is for
# the total number of genes that exist, so we wil cache that many.  If
# the client uses the binomial, then they only ever need to calculate
# the factorial for the number of genes passed into findTerms, but as
# the hypergeometric is the default, we will go with that.
#
# Since :
#
#     n!  = n * (n-1) * (n-2) ... * 1
#
# Then :
#
# log(n!) = log(n * (n-1) * (n-2) ... * 1)
#
#         = log(n) + log(n-1) + log(n-2) ... + log(1)
#

    my ($self) = @_;

    my @logFactorials = (0, 0); # cache of log factorials, initialize for 0 and 1

    my $num = $self->__totalNumGenes;

    for (my $i = 2; $i <= $num; $i++){

	$logFactorials[$i] = $logFactorials[$i-1] + log($i);
	
    }

    # now store the factorials

    $self->{$kLogFactorials} = \@logFactorials;

}



#
# PUBLIC INSTANCE METHODS
#

#####################################################################
sub findTerms{
#####################################################################
# This method returns an array of hash references that indicates what
# terms can annotate the list of genes with what p-value.  The p-value
# is calculated using the hypergeometric distribution (which uses
# sampling without replacement) by default, as the probability of
# seeing that level of annotation to a node or better.  the 'method'
# argument may be used with a value of 'binomial' to use the binomial
# distribution instead, which used sampling with replacement.  This is
# less accurate than the hypergeometric distribution, though for very
# large total numbers in the genome, the differences will be small.
# This option is provided because the binomial distribution is qicker
# to calculate.  
#
# The contents of the hashes in the returned array are:
#
#    key                   value
#    -------------------------------------------------------------------------
#    NODE                  A GO::Node
#
#    PVALUE	           The P-value for having the observed number of
#                          annotations that the provided list of genes
#                          has to that node
#
#    CORRECTED_PVALUE      The CORRECTED_PVALUE is the PVALUE multiplied
#                          by the number of nodes in the minimal set
#                          of hypotheses from which all other
#                          hypotheses can be generated.  A hypothesis
#                          is any node to which 2 or more genes in the
#                          supplied list are annotated, either
#                          directly or indirectly.  The minimal subset
#                          of hypotheses from which all others can be
#                          constructed consists of the union of the
#                          following classes of node: 
#
#                           1).  Leaf hypotheses (i.e. hypotheses which 
#                                have no children that were tested as hypotheses).
#
#                           2).  Hypotheses that have at least one non-hypothesis 
#                                child with an annotation.
#      
#                           3).  Hypotheses with direct annotation.
#
#    NUM_ANNOTATIONS       The number of genes within the provided list that
#                          are annotated to the node.
#
#    TOTAL_NUM_ANNOTATIONS The number of genes across the genome
#                          annotated to the node
#
# The entries are sorted by increasing p-value (ie least likely is first).
#
# It expects to be passed, by reference, a list of gene names for
# which terms will be found.  If a passed in name is ambiguous (see
# documentation for AnnotationProvider), then the following will
# occur:
#
# 1) If the name can be used as a standard name, it will
#    assume that it is that.
#
# 2) Otherwise it will not use it.
#
# If a gene name is not recognized at all, then it will simply be
# ignored, and thus the list will be smaller - I am not yet convinced
# this is the best policy.....
#
# Usage:
#
#    my @pvalueStructures = $termFinder->findTerms(genes=>\@genes);
#
#    my $hypothesis = 1;
#
#    foreach my $pvalue (@pvalueStructures){
#
#    print "-- $hypothesis of ", scalar @pvalueStructures, "--\n",
#
#	"GOID\t", $pvalue->{NODE}->goid, "\n",
#
#	"TERM\t", $pvalue->{NODE}->term, "\n",
#
#	"P-VALUE\t", $pvalue->{PVALUE}, "\n",
#
#	"CORRECTED P-VALUE\t", $pvalue->{CORRECTED_PVALUE}, "\n",
#	
#	"NUM_ANNOTATIONS\t", $pvalue->{NUM_ANNOTATIONS}, " (of ", $pvalue->{TOTAL_NUM_ANNOTATIONS}, ")\n\n";
#
#       $hypothesis++;
#
#    }
#

    my ($self, %args) = @_;

    if (!exists ($args{'genes'})){

	die "You must provide a genes argument";

    }elsif (!defined ($args{'genes'})){

	die "Your genes argument is undefined";

    }

    $self->{$kMethod} = $args{'method'} || 'hypergeometric';

    if (!exists $kAllowedMethods{$self->{$kMethod}}){

	die "$self->{$kMethod} is not an allowed method.  Use one of :". join(", ", keys %kAllowedMethods);

    }

    # what we want to do now, is build up an array of identifiers where
    # that are unambiguous - ie databaseId's
    #
    # This means that when retrieving GOID's, we can always retrieve
    # them by databaseId, which is unambiguous.

    $self->__determineDatabaseIdsFromGenes($args{'genes'});

    if (scalar ($self->__databaseIds) > $self->__totalNumGenes){

	print "You have provided a list corresponding to ", scalar ($self->__databaseIds), "genes, ",
	
	"yet you have indicated that there are only ", $self->__totalNumGenes, " in the genome.\n";

	print "No probabilities can be calculated.\n";

	return (); # simply return an empty list

    } 

    $self->{$kGoCounts} = $self->__buildHashRefOfAnnotations([$self->__databaseIds]);

    $self->__calculatePValues;

    # now what we want to do is calculate pvalues that are corrected for
    # multiple hypothesis testing

    $self->__correctPvalues;

    return $self->__pValues;

}

#
# PRIVATE INSTANCE METHODS
#

#####################################################################
sub __totalNumAnnotatedGenes{
#####################################################################
# This private method returns the number of genes that have any annotation,
# as determined from the AnnotationProvider.  This is set during object
# initialization.

    return $_[0]->{$kTotalNumAnnotatedGenes};

}

#####################################################################
sub __databaseIds{
#####################################################################
# This method returns an array of databaseIds corresponding to the
# genes that were used to initialize the object.

    return @{$_[0]->{$kDatabaseIds}};

}

#####################################################################
sub __pValues{
#####################################################################
# This method returns an array of pValues structures

    return @{$_[0]->{$kPvalues}};

}

#####################################################################
sub __method{
#####################################################################
# This method returns the method by which the client has chosen to
# have their p-values calculated - either binomial or hypergeometric.

    return $_[0]->{$kMethod};

}

#####################################################################
sub __determineDatabaseIdsFromGenes{
#####################################################################
# This method determines a list of databaseIds for the list of
# supplied list of genes for which the client wants to find GO terms.
# It then stores them within the object.
#
# If more than one gene maps to the same databaseId, then the
# databaseId is only put in the list once, and a warning is printed.
#
# If a gene does not map to a databaseId, then an undef is put in the
# list - however, if the same gene name, which does not map to a
# databaseId, is used twice then it will produce only one undef in the
# list.
#
# In addition, it removes leading and trailing whitespace from supplied
# gene names (assuming they should have none) and will skip any names that
# are either empty, or whitespace only.

    my ($self, $genesRef) = @_;

    my (@databaseIds, $databaseId, %databaseIds, %genes);

    foreach my $gene (@{$genesRef}){

	$gene =~ s/^\s+//;
	$gene =~ s/\s+$//;

	next if $gene eq ""; # skip empty names

	if (exists ($genes{$gene})){

	    print "The gene name '$gene' was used more than once.\n";
	    print "It will only be considered once.\n\n";
	    
	    next; # just skip to the next supplied gene

	}

	if ($self->__annotationProvider->nameIsAmbiguous($gene)){

	    print "$gene is an ambiguous name.\n";

	    if ($self->__annotationProvider->nameIsStandardName($gene)){

		print "Since $gene is used as a standard name, it will be assumed to be one.\n\n";
	
		$databaseId = $self->__annotationProvider->databaseIdByStandardName($gene);
	
		push (@databaseIds, $databaseId);
		
	    }else{
		
		print "Since $gene is an ambiguous alias, it will not be used.\n\n";
		
	    }
	    
	}else{

	    # note, if the gene has no annotation, then we will put an
	    # undef into this array of databaseIds - we have to make
	    # sure we deal with this later when getting annotations.	    

	    $databaseId = $self->__annotationProvider->databaseIdByName($gene);

	    # if the total number of genes is equal to the number of
	    # things with some annotation, then there should be no
	    # genes that do not return a databaseId.  If this is the
	    # case, we will warn them.

	    if (!defined $databaseId && $self->__totalNumAnnotatedGenes == $self->__totalNumGenes){

		print "\nThe name '$gene' did not correspond to an entry from the AnnotationProvider.\n";
		print "However, the client has indicated that all genes have annotation.\n";
		print "You should probably check that '$gene' is a real name.\n\n";

	    }

	    push (@databaseIds, $databaseId);

	}

	# if we have a databaseId that we've already seen, we want to
	# make sure we only consider it once.

	if (defined ($databaseId) && exists($databaseIds{$databaseId})){

	    print "More that one gene maps to the same databaseId.\n";
	    print "$gene maps to $databaseId, as did $databaseIds{$databaseId}.\n";
	    print "Only one will be used.\n\n";

	    pop (@databaseIds); # get rid of the extra

	}

	$databaseIds{$databaseId} = $gene if (defined ($databaseId));
	$genes{$gene}             = undef;

    }

    # now store them within the self object

    $self->{$kDatabaseIds} = \@databaseIds;

}

############################################################################
sub __buildHashRefOfAnnotations{
############################################################################
# This private method takes a reference to an array of databaseIds and
# calculates the level of annotations for all GO nodes that those
# databaseIds have either direct or indirect annotation for.  It
# returns a reference to a hash of GO node counts, with the goids
# being the keys, and the number of annotations they have from the
# list of databaseId's being the values.

    my ($self, $databaseIdsRef) = @_;

    my (%goNodeCounts, $goid);

    foreach my $databaseId (@{$databaseIdsRef}) {

	# get goids, if the databaseId is defined

	my @goids = $self->__allGOIDsForDatabaseId($databaseId) if (defined $databaseId);

	if (!@goids) { 
	    
	    # If gene has no annotation, annotate it to the top node
	    # (Gene_Ontology), and its immediate child (the aspect itself)
	    # and the 'unannotated' node.

	    push (@goids, (($self->__ontologyProvider->rootNode->childNodes())[0]->goid, 
			   $self->__ontologyProvider->rootNode->goid,
			   $kUnannotatedNode->goid));
	    
	}

	# increment count for all goids appearing in @goids;

	foreach $goid (@goids) {

	    $goNodeCounts{$goid}++;

	}
    }

    return \%goNodeCounts;

}

############################################################################
sub __allGOIDsForDatabaseId{
############################################################################
# This method returns an array of all GOIDs to which a databaseId is
# annotated, whether explicitly, or implicitly, by virtue of the GO
# node being an ancestor of an explicitly annotated one.  The returned
# array contains no duplicates.
#

    my ($self, $databaseId) = @_;
    
    # generate list of GOIDs if not cached
    
    if (!exists($self->{$kGOIDsForDatabaseIds}->{$databaseId})) {
	
	my %goids; # so we keep the list unique
	
	foreach my $goid (@{$self->__annotationProvider->goIdsByDatabaseId(databaseId => $databaseId,
									   aspect     => $self->__aspect)}) {	    

	    # just in case an annotation is to a goid not present in the ontology

	    if (!$self->__ontologyProvider->nodeFromId($goid)){ # 

		print "$goid, used to annotate $databaseId, does not appear in the ontology.\n";
		
		# don't record any annotations for this databaseId - 
		# will mean it just gets Gene_Ontology, its child, and unannotated

		next;

	    }

	    # record the goid and its ancestors

	    $goids{$goid} = undef;

	    foreach my $ancestor ($self->__ontologyProvider->nodeFromId($goid)->ancestors){

		$goids{$ancestor->goid} = undef;

	    }

	}    
	
	# cache the value
	
	$self->{$kGOIDsForDatabaseIds}->{$databaseId} = [keys %goids];
	
    }
    
    return (@{$self->{$kGOIDsForDatabaseIds}->{$databaseId}});
    
}

#####################################################################
sub __calculatePValues{
#####################################################################
# This method actually determines the p-values of the various levels
# of annotation for the particular GO nodes, and stores them within
# the object.

    my $self = shift;

    my $numDatabaseIds    = scalar $self->__databaseIds;

    my @pvalueArray;
    
    foreach my $goid ($self->__allGoIdsForList) {
	
	next if ($self->__numAnnotationsToGoId($goid) == 1); # skip GO nodes with only one gene annotation

	push (@pvalueArray, $self->__processOneGOID($goid, $numDatabaseIds));

    }

    @pvalueArray = sort {$a->{PVALUE} <=> $b->{PVALUE}} @pvalueArray;

    $self->{$kPvalues} = \@pvalueArray

}

############################################################################
sub __processOneGOID{
############################################################################
# This processes one GOID.  It determines the number of annotations to
# the current GOID, and the P-value of that number of annotations.
# It returns a hash reference encoding that information.

    my ($self, $goid, $numGenes) = @_;

    my $totalNumAnnotationsToGoId = $self->__totalNumAnnotationsToGoId($goid);
    my $numAnnotationsToGoId      = $self->__numAnnotationsToGoId($goid);
    my $totalNumGenes             = $self->__totalNumGenes();
    my $p                         = $totalNumAnnotationsToGoId / $totalNumGenes;
    my $method                    = $self->__method;

    my $pvalue = 0;
    
    for (my $j = $self->__numAnnotationsToGoId($goid); $j <= $numGenes; $j++) {
	
	if ($method eq 'hypergeometric'){

	    $pvalue += $self->__hypergeometric($j, $numGenes, $totalNumAnnotationsToGoId, $totalNumGenes);

	}else{

	    $pvalue += $self->__binomial($j, $numGenes, $p);

	}

    }
    
    my $node = $self->__ontologyProvider->nodeFromId($goid) || $kUnannotatedNode;
    
    my $hashRef = {
	
	NODE                  => $node,
	PVALUE		      => $pvalue,
	NUM_ANNOTATIONS       => $numAnnotationsToGoId,
	TOTAL_NUM_ANNOTATIONS => $totalNumAnnotationsToGoId
	    
	};
    
    return $hashRef;

}

############################################################################
sub __numAnnotationsToGoId{
############################################################################
# This private method returns the number of annotations to a
# particular GOID for the list of genes supplied to the findTerms
# method.

    my ($self, $goid) = @_;

    return $self->{$kGoCounts}->{$goid};

}

############################################################################
sub __totalNumAnnotationsToGoId{
############################################################################
# This returns the total number of genes that have been annotated to a
# particular GOID based on all annotations.

    my ($self, $goid) = @_;

    return $self->{$kTotalGoNodeCounts}->{$goid};
}

############################################################################
sub __totalNumGenes{
############################################################################
# This returns the total number of genes that exist for the particular 
# organism in question. Unannotated genes are included in this count.


    return $_[0]->{$kArgs}{totalNumGenes};

}

############################################################################
sub __binomial{
############################################################################
# This method calculates the binomial distribution probability, given
# $j successes from $numGenes trials, given a probability of $p of
# success of any given trial.  Trials are assumed to be independent.
#
#
# The binomial distribution probability is equal to:
#
#    ( numGenes choose j ) * p^j * (1-p)^(numGenes - j)
#
#
# we can do this in log space, to avoid buffer overflows:
#
#   log(( numGenes choose j ) * p^j * (1-p)^(numGenes - j))
#
# = log(numGenes choose j) + log(p^j) + log((1-p)^(numGenes - j))
#
# = log(numGenes choose j) + j*log(p) + (numGenes - j) * log(1-p)
#
# Note, we have a problem with log of zero if p is 1.  In this case,
# the binomial is equal to 1, because if the probability of a success
# is 1 then all trials will be successful, so the probability of j of
# numGenes successes must be 1.
#

    my ($self, $j, $numGenes, $p) = @_;

    if ($p != 1){

	return exp($self->__logNCr($numGenes, $j) + $j * log($p) + ($numGenes - $j) * log(1-$p));

    }else{

	return 1;

    }	

}

############################################################################
sub __hypergeometric{
############################################################################
# This method returns the hypergeometric probability value for
# sampling without replacement.  The calculation is the probability of
# picking x positives from a sample of n, given that there are M
# positives in a population of N.
#
# The value is calculated as:
#
#       (M choose x) (N-M choose n-x)
# P =   -----------------------------
#               N choose n
#
# where generically n choose r is number of permutations by which r
# things can be chosen from a population of n (see __logNCr())
#
# However, given that these n choose r values may be extremely high (as they are
# are calculated using factorials) it is safer to do this instead in log space,
# as we are far less likely to have an overflow.
#
# thus :
#
# log(P) = log(M choose x) + log(N-M choose n-x) - log (N choose n);
#
# this means we can now calculate log(n choose R) for our
# hypergeometric calculation (see below).
#
    my ($self, $x, $n, $M, $N) = @_;

    return exp($self->__logNCr($M, $x) + $self->__logNCr($N - $M, $n-$x) - $self->__logNCr($N, $n));

}

############################################################################
sub __logNCr{
############################################################################
# This method returns the log of n choose R.  This means that it can do the
# calculation in log space itself.
#
#
#           n!
# nCr =  ---------
#        r! (n-r)!
#
# which means:
#
#
#
# log(nCr) = log(n!) - (log(r!) + log((n-r)!))
#

    my ($self, $n, $r) = @_;

    if (!exists $self->{$kLogNCr}{$n}{$r}){

	$self->{$kLogNCr}{$n}{$r} = $self->__logFact($n) - ($self->__logFact($r) + $self->__logFact($n - $r));

    }

    return $self->{$kLogNCr}{$n}{$r};

}

############################################################################
sub __logFact{
############################################################################
# This method returns the log of a factorial, from our previously calculated
# cache (see __cacheLogFactorials()).
#

    return $_[0]->{$kLogFactorials}->[$_[1]];

}    

############################################################################
sub __allGoIdsForList{
############################################################################
# This returns an array of GOIDs to which genes in the passed in gene
# list were directly or indirectly annotated.

    return keys %{$_[0]->{$kGoCounts}};

}

############################################################################
sub __correctPvalues{
############################################################################
# This method corrects the pvalues for multiple hypothesis testing.
# Because each hypothesis (a GO node) is not an independent
# hypothesis, then full correction (multiplying by the number of nodes
# (hypotheses) considered) is not appropriate.  Instead, correction is
# done by multiplying each pvalue by the number of nodes in the
# minimal subset from which all hypotheses can be reconstructed.  This
# is the number of nodes whose level of annotatation cannot be solely
# reconstructed from child nodes that were hypotheses themselves.
# This corresponds to the union of the following three classes of node:
#
# 1).  Leaf hypotheses (i.e. hypotheses which have no children that were tested as
#      hypotheses).
#
# 2).  Hypotheses that have at least one non-hypothesis child with an annotation.
#      
# 3).  Hypotheses with direct annotation.

    my ($self) = @_;

    my @leafHypotheses = $self->__leafHypotheses;

    my @hypothesesWithNonHypothesisAnnotatedChildren = $self->__hypothesesWithNonHypothesisAnnotatedChildren;

    my @directlyAnnotatedHypotheses = $self->__directlyAnnotatedHypotheses;

    # now determine the number of unique hypotheses that are needed to
    # reconstruct all hypotheses

    my %goids;

    foreach my $goid (@leafHypotheses, @hypothesesWithNonHypothesisAnnotatedChildren, @directlyAnnotatedHypotheses){
	
	$goids{$goid} = undef;
	
    }

    my $minimalNumberOfHypotheses = scalar keys %goids;

    # now correct pvalues, by multiplying by $minimalNumberOfHypotheses
    
    foreach my $hypothesis ($self->__pValues){

	$hypothesis->{CORRECTED_PVALUE} = $hypothesis->{PVALUE} * $minimalNumberOfHypotheses;

	# make sure we have a ceiling of 1

	$hypothesis->{CORRECTED_PVALUE} = 1 if ($hypothesis->{CORRECTED_PVALUE} > 1);

    }

}

############################################################################
sub __leafHypotheses{
############################################################################
# This method returns the goids of hypothesis nodes that don't
# have children that were hypothesis nodes.  For each hypothesis node
# we just have to check that none of its children had two or more
# annotations, as that was the requirement to be considered as
# hypothesis.

    my ($self) = @_;

    my @leafHypotheses;

  PARENT:

    foreach my $hypothesis ($self->__pValues){
	
	# now go through each child of this goid, and see if it had 2
	# or more annotations
	
	foreach my $childNode ($hypothesis->{NODE}->childNodes){
	    
	    # if any of the children have more than one annotation, we skip
	    # this parent entirely

	    next PARENT if (defined($self->__numAnnotationsToGoId($childNode->goid)) &&

			    $self->__numAnnotationsToGoId($childNode->goid) > 1);

	}

	# if we get here, either the considered node had no children,
	# or none of its children were hypotheses (had 2 or more annotations from our list)

	push (@leafHypotheses, $hypothesis->{NODE}->goid);

    }

    return @leafHypotheses;

}

############################################################################
sub __hypothesesWithNonHypothesisAnnotatedChildren{
############################################################################
# This method returns an array of goids that correspond to hypotheses that
# were tested, that have child nodes with annotation that were not tested
# as hypotheses.  Such children would only have a single annotation.

    my ($self) = @_;

    my @hypotheses;

    foreach my $hypothesis ($self->__pValues){
	
	# now go through each child of this goid, and see if it had
	# only 1 annotation
	
	foreach my $childNode ($hypothesis->{NODE}->childNodes){
	    
	    # if any of the children have exactly one annotation
	    # we record the node, and move on to looking at the next
	    # hypothesis

	    if (defined($self->__numAnnotationsToGoId($childNode->goid)) &&
		
		$self->__numAnnotationsToGoId($childNode->goid) == 1){

		push (@hypotheses, $hypothesis->{NODE}->goid);

		last; # don't need to check the rest of the children

	    }

	}

    }

    return @hypotheses;

}

############################################################################
sub __directlyAnnotatedHypotheses{
############################################################################
# This method returns an array of all the hypotheses that are directly
# annotated themselves.

    my ($self) = @_;

    my %directlyAnnotatedNodes;

    # first we have to work out which nodes are directly annotated by
    # our list of genes.  We could have done this earlier, but rather
    # than muddy the __buildHashRefOfAnnotations and __allGOIDsForDatabaseId
    # code, we'll do it more cleanly here

    foreach my $databaseId ($self->__databaseIds) {

	next if !defined $databaseId; # skip those with no databaseId

	foreach my $goid (@{$self->__annotationProvider->goIdsByDatabaseId(databaseId => $databaseId,
									   aspect     => $self->__aspect)}) {

	    $directlyAnnotatedNodes{$goid} = undef;

	}

    }

    # now check which of out hypotheses are directly annotated

    my @directlyAnnotatedHypotheses;

    foreach my $hypothesis ($self->__pValues){

	# if the hypothesis is directly annotated

	if (exists ($directlyAnnotatedNodes{$hypothesis->{NODE}->goid})){

	    # record it

	    push (@directlyAnnotatedHypotheses, $hypothesis->{NODE}->goid);

	}

    }

    return (@directlyAnnotatedHypotheses);
    
}

############################################################################
sub __annotationProvider{
############################################################################
# This private method returns the annotationProvider that was used
# during construction.

    return $_[0]->{$kArgs}{annotationProvider};

}

############################################################################
sub __ontologyProvider{
############################################################################
# This private methid returns the ontologyProvider that was used
# during construction.

    return $_[0]->{$kArgs}{ontologyProvider};

}

############################################################################
sub __aspect{
############################################################################

    return $_[0]->{$kArgs}{aspect};

}

1; # to make perl happy


__END__

#####################################################################
#
#  POD Documentation from here on down
#
#####################################################################

=pod

=head1 Instance Constructor

=head2 new

This is the constructor.  It expects to be passed named arguments for
an annotationProvider, an ontologyProvider, and how many genes in
total exist.  In addition, it must be told the aspect of the ontology
provider, so that it knows how to query the annotationProvider.

Usage :

    my $termFinder = GO::TermFinder->new(annotationProvider=> $annotationProvider,
					 ontologyProvider  => $ontologyProvider,
					 totalNumGenes     => $num,
					 aspect            => <P|C|F>);

=head1 Instance Methods

=head2 findTerms

This method returns an array of hash references that indicates what
terms can annotate the list of genes with what P-value.  The
contents of the hashes in the returned array are:

    key                   value
    -------------------------------------------------------------------------
    NODE                  A GO::Node

    PVALUE		  The P-value for having the observed number of
                          annotations that the provided list of genes
                          has to that node

    CORRECTED_PVALUE      The CORRECTED_PVALUE is the PVALUE multiplied
                          by the number of nodes in the minimal set
                          of hypotheses from which all other
                          hypotheses can be generated.  A hypothesis
                          is any node to which 2 or more genes in the
                          supplied list are annotated, either
                          directly or indirectly.  The minimal subset
                          of hypotheses from which all others can be
                          constructed consists of the union of the
                          following classes of node: 

                           1).  Leaf hypotheses (i.e. hypotheses which 
                                have no children that were tested as hypotheses).

                           2).  Hypotheses that have at least one non-hypothesis 
                                child with an annotation.
      
                           3).  Hypotheses with direct annotation.

    NUM_ANNOTATIONS       The number of genes within the provided list that
                          are annotated to the node.

    TOTAL_NUM_ANNOTATIONS The number of genes across the genome
                          annotated to the node

The entries are sorted by increasing p-value (ie least likely is
first).  

This method expects to be passed, by reference, a list of gene names
for which terms will be found.  If a passed in name is ambiguous (see
AnnotationProvider), then the following will occur:


    1) If the name can be used as a standard name, it will assume that
       it is that.

    2) Otherwise it will not use it.

Currently a warning will be printed to STDOUT in the case of an
ambiguous name being used.

The passed in gene names are converted into a list of databaseIds.  If
a gene does not map to a databaseId, then an undef is put in the list
- however, if the same gene name, which does not map to a databaseId,
is used twice then it will produce only one undef in the list.  If
more than one gene name maps to the same databaseId (either because
you used the same name twice, or you used an alias as well), then that
databaseId is only put into the list once, and a warning is printed.


If a gene name does not have any information returned from the
AnnotationProvider, then it is assumed that the gene is entirely
unannotated.  For these purposes, TermFinder annotates such genes to
the root node (Gene_Ontology), its immediate child (which indicates
the aspect of the ontology (such as biological_process), and a dummy
go node, corresponding to unannotated.  This node will have a goid of
'GO:XXXXXXX', and a term name of 'unannotated'.  No other information
will be set up for this GO::Node, so you should not count on being
able to retrieve it.  What it does mean is that you can determine if
the predominant feature of a set of genes is that they have no
annotation.

If more genes are provided that have been indicated exist in the
genome (as provided during object construction), then an error message
will be printed out, and an empty list will be returned.

Usage:

    my @pvalueStructures = $termFinder->findTerms(genes=>\@genes);

    my $hypothesis = 1;						    

    foreach my $pvalue (@pvalueStructures){

    print "-- $hypothesis of ", scalar @pvalueStructures, "--\n",

	"GOID\t", $pvalue->{NODE}->goid, "\n",

	"TERM\t", $pvalue->{NODE}->term, "\n",

	"P-VALUE\t", $pvalue->{PVALUE}, "\n",

	"CORRECTED P-VALUE\t", $pvalue->{CORRECTED_PVALUE}, "\n",
	
	"NUM_ANNOTATIONS\t", $pvalue->{NUM_ANNOTATIONS}, " (of ", $pvalue->{TOTAL_NUM_ANNOTATIONS}, ")\n\n";

    $hypothesis++;

    }

=head1 Authors

    Gavin Sherlock; sherlock@genome.stanford.edu
    Elizabeth Boyle; ell@mit.edu

=cut
