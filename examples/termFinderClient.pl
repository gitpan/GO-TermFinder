#!/usr/bin/perl

# $Id: termFinderClient.pl,v 1.5 2003/12/03 02:33:19 sherlock Exp $

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
    
    "NUM_ANNOTATIONS\t", $pvalue->{NUM_ANNOTATIONS}, " (of ", $pvalue->{TOTAL_NUM_ANNOTATIONS}, ")\n",

    "ANNOTATED_GENES\t", join(", ", values (%{$pvalue->{ANNOTATED_GENES}})), "\n\n";
    
    $hypothesis++;
    
}

=pod

=head1 NAME

termFinderClient.pl - interactive client to find significant GO terms for a list of genes

=head1 SYNOPSIS

This program is a very simply client for the GO::TermFinder object,
that prompts a user for the various pieces of information that are
required to determine significant GO terms associated with the list of
genes.  It uses a p-value cut-off (for the corrected p-value) of .001,
and simply prints the information back to the screen.  An example
below uses the genes : YPL250C, MET11, MXR1, MET17, SAM3, MET28, STR3,
MMP1, MET1, YIL074C, MHT1, MET14, MET16, MET3, MET10, ECM17, MET2,
MUP1 and MET6, which formed the so called methionine cluster from
Spellman et al, 1998:

    > termFinderClient.pl
    Enter the fully qualified name of your ontology file:
    process.ontology
    What is the aspect of this ontology (F, P, or C)?
    P
    Enter the fully qualified name of your associations file:
    gene_association.sgd
    Enter a the fully qualified name of your file with a list of genes for which to find term:
    genes.txt
    How many genes (roughly) exist within the organism?
    7300
    Finding terms...
    -- 1 of 28--
    GOID    GO:0006790
    TERM    sulfur metabolism
    P-VALUE 4.72498557073697e-25
    CORRECTED P-VALUE       6.14248124195806e-24
    NUM_ANNOTATIONS 13 (of 51)
    
    -- 2 of 28--
    GOID    GO:0000096
    TERM    sulfur amino acid metabolism
    P-VALUE 3.29201306923395e-23
    CORRECTED P-VALUE       4.27961699000414e-22
    NUM_ANNOTATIONS 11 (of 29)
    
    -- 3 of 28--
    GOID    GO:0006555
    TERM    methionine metabolism
    P-VALUE 9.48135258251747e-20
    CORRECTED P-VALUE       1.23257583572727e-18
    NUM_ANNOTATIONS 9 (of 20)
    
    -- 4 of 28--
    GOID    GO:0006520
    TERM    amino acid metabolism
    P-VALUE 1.06665546159811e-16
    CORRECTED P-VALUE       1.38665210007754e-15
    NUM_ANNOTATIONS 12 (of 145)
    
      -- 5 of 28--
    GOID    GO:0009066
    TERM    aspartate family amino acid metabolism
    P-VALUE 2.44942452266491e-16
    CORRECTED P-VALUE       3.18425187946438e-15
    NUM_ANNOTATIONS 9 (of 42)
    
    -- 6 of 28--
    GOID    GO:0006519
    TERM    amino acid and derivative metabolism
    P-VALUE 2.84201383187083e-16
    CORRECTED P-VALUE       3.69461798143208e-15
    NUM_ANNOTATIONS 12 (of 157)
    
    -- 7 of 28--
    GOID    GO:0009308
    TERM    amine metabolism
    P-VALUE 6.06713484278472e-16
    CORRECTED P-VALUE       7.88727529562014e-15
    NUM_ANNOTATIONS 12 (of 167)
    
    -- 8 of 28--
    GOID    GO:0000103
    TERM    sulfate assimilation
    P-VALUE 3.61076971610707e-15
    CORRECTED P-VALUE       4.69400063093919e-14
    NUM_ANNOTATIONS 6 (of 8)
    
    -- 9 of 28--
    GOID    GO:0006791
    TERM    sulfur utilization
    P-VALUE 3.61076971610707e-15
    CORRECTED P-VALUE       4.69400063093919e-14
    NUM_ANNOTATIONS 6 (of 8)
    
    -- 10 of 28--
    GOID    GO:0006082
    TERM    organic acid metabolism
    P-VALUE 2.41612157270023e-14
    CORRECTED P-VALUE       3.1409580445103e-13
    NUM_ANNOTATIONS 12 (of 226)
    
    -- 11 of 28--
    GOID    GO:0019752
    TERM    carboxylic acid metabolism
    P-VALUE 2.41612157270023e-14
    CORRECTED P-VALUE       3.1409580445103e-13
    NUM_ANNOTATIONS 12 (of 226)
    
    -- 12 of 28--
    GOID    GO:0000097
    TERM    sulfur amino acid biosynthesis
    P-VALUE 3.75638037739466e-12
    CORRECTED P-VALUE       4.88329449061306e-11
    NUM_ANNOTATIONS 5 (of 8)
    
    -- 13 of 28--
    GOID    GO:0008652
    TERM    amino acid biosynthesis
    P-VALUE 1.25659010245536e-07
    CORRECTED P-VALUE       1.63356713319197e-06
    NUM_ANNOTATIONS 6 (of 99)
    
    -- 14 of 28--
    GOID    GO:0000101
    TERM    sulfur amino acid transport
    P-VALUE 1.49023707732037e-07
    CORRECTED P-VALUE       1.93730820051648e-06
    NUM_ANNOTATIONS 3 (of 5)
    
    -- 15 of 28--
    GOID    GO:0009086
    TERM    methionine biosynthesis
    P-VALUE 1.49023707732037e-07
    CORRECTED P-VALUE       1.93730820051648e-06
    NUM_ANNOTATIONS 3 (of 5)
    
    -- 16 of 28--
    GOID    GO:0009309
    TERM    amine biosynthesis
    P-VALUE 1.89250974784783e-07
    CORRECTED P-VALUE       2.46026267220218e-06
    NUM_ANNOTATIONS 6 (of 106)
    
    -- 17 of 28--
    GOID    GO:0009067
    TERM    aspartate family amino acid biosynthesis
    P-VALUE 1.19029139172751e-05
    CORRECTED P-VALUE       0.000154737880924577
    NUM_ANNOTATIONS 3 (of 18)
    
    -- 18 of 28--
    GOID    GO:0006865
    TERM    amino acid transport
    P-VALUE 5.80646562760316e-05
    CORRECTED P-VALUE       0.000754840531588411
    NUM_ANNOTATIONS 3 (of 30)

=head1 AUTHORS

Gavin Sherlock, sherlock@genome.stanford.edu

=cut
