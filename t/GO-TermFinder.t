use Test;
BEGIN { plan tests => 27 };

# File       : GO-TermFinder.t
# Author     : Gavin Sherlock
# Date Begun : September 1st 2003

# $Id: GO-TermFinder.t,v 1.1 2003/10/16 17:23:35 sherlock Exp $

# This file forms a set of tests for the GO::TermFinder class

use strict;
use warnings;
use diagnostics;

use GO::TermFinder;
use GO::AnnotationProvider::AnnotationParser;
use GO::OntologyProvider::OntologyParser;

my $ontologyFile   = "t/process.ontology";
my $annotationFile = "t/gene_association.sgd"; 

# we'll check that all the public methods still exist in the interface

my @methods = qw (findTerms);

my $ontology = GO::OntologyProvider::OntologyParser->new(ontologyFile=>$ontologyFile);

my $annotation = GO::AnnotationProvider::AnnotationParser->new(annotationFile=>$annotationFile);

my $termFinder = GO::TermFinder->new(annotationProvider=> $annotation,
				     ontologyProvider  => $ontology,
				     totalNumGenes     => 6905,
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

# the following are what should be the 10 most significant goids for
# this set of genes using this frozen dataset

my @topGoids = ("GO:0006790",
		"GO:0000096",
		"GO:0006555",
		"GO:0009066",
		"GO:0000103",
		"GO:0006791",
		"GO:0006520",
		"GO:0006519",
		"GO:0009308",
		"GO:0019752");

# now check that these are returned by the TermFinder

for (my $i = 0; $i< @topGoids; $i++){

    ok($pvalues[$i]->{NODE}->goid, $topGoids[$i]);

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
