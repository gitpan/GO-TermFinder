#!/usr/bin/perl

# $Id: termFinderClient.pl,v 1.2 2003/03/03 18:25:50 sherlock Exp $

# Date   : 3rd February 2003
# Author : Gavin Sherlock

# This program is a very simply client for the GO::TermFinder object,
# that prompts a user for the various pieces of information that are
# required to generate the terms, and simply prints the information
# back to the screen.

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

use strict;
use warnings;
use diagnostics;

use lib "../lib";
use GO::TermFinder;
use GO::AnnotationProvider::AnnotationParser;
use GO::OntologyProvider::OntologyParser;

print "Enter the fully qualified name of your ontology file:\n";

chomp(my $ontologyFile = <STDIN>);

print "What is the aspect of this ontology (F, P, or C)?\n";

chomp(my $aspect = uc(<STDIN>));

print "Enter the fully qualified name of your associations file:\n";

chomp(my $annotationFile = <STDIN>); 

print "Enter a the fully qualified name of your file with a list of genes for which to find term:\n";

chomp(my $genesFile = <STDIN>);

print "How many genes (roughly) exist within the organism?\n";

chomp(my $numGenes = <STDIN>);

print "Finding terms...\n";

my @genes;

open (IN, $genesFile) || die "Cannot open $genesFile :$!";

while (<IN>){

    chomp;

    push (@genes, $_);

}

close IN;

my $ontology   = GO::OntologyProvider::OntologyParser->new(ontologyFile=>$ontologyFile);

my $annotation = GO::AnnotationProvider::AnnotationParser->new(annotationFile=>$annotationFile);

my $termFinder = GO::TermFinder->new(annotationProvider=> $annotation,
				     ontologyProvider  => $ontology,
				     totalNumGenes     => $numGenes,
				     aspect            => $aspect);

my @pvalues    = $termFinder->findTerms(genes=>\@genes);

# now just print the info back to the client

my $hypothesis = 1;

foreach my $pvalue (@pvalues){
    
    next if ($pvalue->{CORRECTED_PVALUE} > 0.001);
    
    print "-- $hypothesis of ", scalar @pvalues, "--\n",
    
    "GOID\t", $pvalue->{NODE}->goid, "\n",
    
    "TERM\t", $pvalue->{NODE}->term, "\n",
    
    "P-VALUE\t", $pvalue->{PVALUE}, "\n",
    
    "CORRECTED P-VALUE\t", $pvalue->{CORRECTED_PVALUE}, "\n",
    
    "NUM_ANNOTATIONS\t", $pvalue->{NUM_ANNOTATIONS}, " (of ", $pvalue->{TOTAL_NUM_ANNOTATIONS}, ")\n\n";
    
    $hypothesis++;
    
}
