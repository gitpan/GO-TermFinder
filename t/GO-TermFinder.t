use Test;
BEGIN { plan tests => 1123 };

# File       : GO-TermFinder.t
# Author     : Gavin Sherlock
# Date Begun : September 1st 2003

# $Id: GO-TermFinder.t,v 1.4 2003/12/03 02:30:25 sherlock Exp $

# This file forms a set of tests for the GO::TermFinder class

use strict;
use warnings;
use diagnostics;

use GO::TermFinder;
use GO::AnnotationProvider::AnnotationParser;
use GO::OntologyProvider::OntologyParser;

$|=1;

my $ontologyFile   = "t/process.ontology";
my $annotationFile = "t/gene_association.sgd"; 

# we'll check that all the public methods still exist in the interface

my @methods = qw (findTerms);

my $ontology = GO::OntologyProvider::OntologyParser->new(ontologyFile=>$ontologyFile);

my $annotation = GO::AnnotationProvider::AnnotationParser->new(annotationFile=>$annotationFile);

my $termFinder = GO::TermFinder->new(annotationProvider=> $annotation,
				     ontologyProvider  => $ontology,
				     aspect            => 'P');

ok($termFinder->isa("GO::TermFinder"));

# check the object returns a code reference when asked it it can do a
# method that should exist

foreach my $method (@methods){

    ok(ref($termFinder->can($method)), "CODE");

}

# now check that the findTerms method actually returns the correct
# answers for a selected list of genes (actually the methionine
# cluster from Spellman et al, 1998).

my @pvalues = $termFinder->findTerms(genes=>[qw(YPL250C
						MET11
						MXR1
						MET17
						SAM3
						MET28
						STR3
						MMP1
						MET1
						YIL074C
						MHT1
						MET14
						MET16
						MET3
						MET10
						ECM17
						MET2
						MUP1
						MET6)]);

&testHypotheses(@pvalues);

# now let's run exactly the same test again, but with a different
# casing of the genes

my @newpvalues = $termFinder->findTerms(genes=>[qw(ypl250c
						   Met11
						   mxr1
						   Met17
						   SAM3
						   met28
						   Str3
						   MMp1
						   mET1
						   YIl074c
						   Mht1
						   mEt14
						   Met16
						   Met3
						   mET10
						   ecm17
						   Met2
						   MuP1
						   MeT6)]);

# and compare that the stuff returned looks exactly the same

&compareHypotheses(\@pvalues, \@newpvalues);

# now let's test the functionality of using a defined population
# to create the background distribution.  If we simply say that
# the defined population is every gene from the annotation parser
# then we should get the same result
    
my $newTermFinder = GO::TermFinder->new(annotationProvider=> $annotation,
					ontologyProvider  => $ontology,
					population        => [$annotation->allDatabaseIds],
					aspect            => 'P');

my @poppvalues = $newTermFinder->findTerms(genes=>[qw(ypl250c
						      Met11
						      mxr1
						      Met17
						      SAM3
						      met28
						      Str3
						      MMp1
						      mET1
						      YIl074c
						      Mht1
						      mEt14
						      Met16
						      Met3
						      mET10
						      ecm17
						      Met2
						      MuP1
						      MeT6)]);

# again, check that the stuff returned looks exactly the same

&compareHypotheses(\@pvalues, \@poppvalues);

# now try using a TermFinder with a limited population of just a few genes.
# All of the returned nodes should have a probability of 1

my $nextTermFinder = GO::TermFinder->new(annotationProvider=> $annotation,
					 ontologyProvider  => $ontology,
					 population        => [qw(ypl250c
								  Met11
								  mxr1
								  Met17
								  SAM3
								  met28
								  Str3
								  MMp1
								  mET1
								  YIl074c
								  Mht1
								  mEt14
								  Met16
								  Met3
								  mET10
								  ecm17
								  Met2
								  MuP1
								  MeT6)],
					 aspect            => 'P');


@pvalues = $nextTermFinder->findTerms(genes=>[qw(ypl250c
						 Met11
						 mxr1
						 Met17
						 SAM3
						 met28
						 Str3
						 MMp1
						 mET1
						 YIl074c
						 Mht1
						 mEt14
						 Met16
						 Met3
						 mET10
						 ecm17
						 Met2
						 MuP1
						 MeT6)]);

foreach my $pvalue (@pvalues){

    # round the pvalue to 2 decimal places, to avoid precision problems

    my $val = sprintf("%.2f", $pvalue->{PVALUE});

    ok($val, "1.00");

}

# Now let's test some of the math function insider the TermFinder.

# These are private functions, but it's important that we test that
# they always return the expected values.

# check that the logfactorial is okay

# first calculate factorials from 0 through 10

my @factorials = (1, 1); # initialize for 0 and 1

for (my $i = 2; $i <= 10; $i++){

    $factorials[$i] = $factorials[$i-1] * $i;

}

# now check them against the log values from TermFinder

for (my $i = 0; $i <= 10; $i++){

    ok(log($factorials[$i]), $termFinder->__logFact($i));

}

# now check that the __logNCr method is working correctly

# test that we get the correct value as if 6 had been chosen out of 10,
# given that:
#
#           n!
# nCr =  ---------
#        r! (n-r)!

{ # lexically scope, to prevent collision between this $n and the $n
  # used below

    my $n = 10;
    my $r = 6;
    
    my $nChooseR = $factorials[$n] / ($factorials[$r] * $factorials[$n-$r]);
    
    # now check against the log value that the TermFinder will return
    
    ok($termFinder->__logNCr($n, $r), log($nChooseR));

}

# now let's test that the hypergeometric function works correctly
#
# we'll do a simple test for the probability of picking 3 out of 5,
# given that in the population there is 4 out of 10
#
# The calculation is the probability of picking x positives from a
# sample of n, given that there are M positives in a population of N.
#
# The value is calculated as:
#
#       (M choose x) (N-M choose n-x)
# P =   -----------------------------
#               N choose n
#

my $M = 4;
my $N = 10;

my $n = 5;
my $x = 3;

my $a = $factorials[$M] / ($factorials[$x] * $factorials[$M-$x]);
my $b = $factorials[$N - $M] / ($factorials[$n - $x] * $factorials[($N - $M) - ($n - $x)]);
my $c = $factorials[$N] / ($factorials[$n] * $factorials[$N-$n]);

my $probability = ($a * $b) / $c;

ok($probability, $termFinder->__hypergeometric($x, $n, $M, $N));

# now we want to check the pvalue using the hypergeometric
#
# the pvalue is the probability of getting x or more from a sample of
# n, given M positives in a population of N
#
# We'll use the same example as above, and calculate the pvalue for 3
# of 5, given 4 of 10 in the population

my $pvalue = 0;

for (my $i = $x; $i <= $n; $i++){

    my $a = $factorials[$M] / ($factorials[$i] * $factorials[$M-$i]);
    my $b = $factorials[$N - $M] / ($factorials[$n - $i] * $factorials[($N - $M) - ($n - $i)]);
    my $c = $factorials[$N] / ($factorials[$n] * $factorials[$N-$n]);

    my $probability = ($a * $b) / $c;

    $pvalue += $probability;

}

# because the TermFinder module uses log space internally to calculate
# factorials and nChooseR, it is not as precise as this test-suite.  Thus
# we need to reduce the precision a little.
#
# Should get TermFinder to use BigInt sometime....

$pvalue = sprintf("%.8f", $pvalue);

my $test = sprintf("%.8f", $termFinder->__pValueByHypergeometric($x, $n, $M, $N));

ok($pvalue, $test);

$pvalue = 0;

for (my $i = 0; $i < $x; $i++){

    my $a = $factorials[$M] / ($factorials[$i] * $factorials[$M-$i]);
    my $b = $factorials[$N - $M] / ($factorials[$n - $i] * $factorials[($N - $M) - ($n - $i)]);
    my $c = $factorials[$N] / ($factorials[$n] * $factorials[$N-$n]);

    my $probability = ($a * $b) / $c;

    $pvalue += $probability;

}

$pvalue = 1 - $pvalue;

$pvalue = sprintf("%.8f", $pvalue);

ok($pvalue, $test);


######################################################################################
sub testHypotheses{
######################################################################################
# the following are what should be the 11 most significant goids for
# the test set of genes using this frozen dataset.  Note if two nodes
# result in the same p-value, they will be sorted by GOID, using a
# text sort.  There are two cases below where this happens.

    my @pvalues = @_;

    my @topGoids = ("GO:0006790",
		    "GO:0000096",
		    "GO:0006555",
		    "GO:0009066",
		    "GO:0000103",  # this and GO:0006791 actually give the same p-value
		    "GO:0006791",
		    "GO:0006520",
		    "GO:0006519",
		    "GO:0009308",
		    "GO:0006082",  # this and GO:0019752 actually give the same p-value
		    "GO:0019752");

    # now check that these are returned by the TermFinder

    for (my $i = 0; $i< @topGoids; $i++){

	ok($pvalues[$i]->{NODE}->goid, $topGoids[$i]);        

    }

}

######################################################################################
sub compareHypotheses{
######################################################################################
# This subroutine expects to receive two arrays (by reference) of hypotheses
# generated by GO::TermFinder.  It will check whether they are identical

    my ($ref1, $ref2) = @_;

    for (my $i = 0; $i < @{$ref1}; $i++){

	ok($ref1->[$i]->{PVALUE},                $ref2->[$i]->{PVALUE});  
	ok($ref1->[$i]->{CORRECTED_PVALUE},      $ref2->[$i]->{CORRECTED_PVALUE});  
	ok($ref1->[$i]->{NUM_ANNOTATIONS},       $ref2->[$i]->{NUM_ANNOTATIONS});  
	ok($ref1->[$i]->{TOTAL_NUM_ANNOTATIONS}, $ref2->[$i]->{TOTAL_NUM_ANNOTATIONS}); 
	ok($ref1->[$i]->{NODE}->goid,            $ref2->[$i]->{NODE}->goid);  
	ok($ref1->[$i]->{NODE}->term,            $ref2->[$i]->{NODE}->term);  

	# now check the genes

	# same number

	ok(scalar keys (%{$ref1->[$i]->{ANNOTATED_GENES}}), scalar keys (%{$ref2->[$i]->{ANNOTATED_GENES}}));

	foreach my $gene (keys (%{$ref1->[$i]->{ANNOTATED_GENES}})){

	    # each one exists

	    ok (exists $ref2->[$i]->{ANNOTATED_GENES}{$gene});

	    # and has the same name - not has to be done case-insensitively
	    
	    ok(uc($ref1->[$i]->{ANNOTATED_GENES}{$gene}), uc($ref2->[$i]->{ANNOTATED_GENES}{$gene}));

	}

    }

}


=head1 Modifications

 List them here.

 CVS information:

 # $Author: sherlock $
 # $Date: 2003/12/03 02:30:25 $
 # $Log: GO-TermFinder.t,v $
 # Revision 1.4  2003/12/03 02:30:25  sherlock
 # added in a bunch of tests to be more precise in the testing of the
 # term finder, and to test the functionality of providing a population
 # of genes from which to calculate the background distribution
 #
 # Revision 1.3  2003/11/22 00:09:12  sherlock
 # added some new tests to see that differently cased versions of the
 # gene names still give the same result.
 #

=cut
