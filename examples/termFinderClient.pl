#!/usr/bin/perl

# $Id: termFinderClient.pl,v 1.6 2004/05/06 01:36:49 sherlock Exp $

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

use GO::TermFinderReport::Text;

use GO::Utils::File qw (GenesFromFile);

print "Enter the fully qualified name of your ontology file:\n";

chomp(my $ontologyFile = <STDIN>);

print "What is the aspect of this ontology (F, P, or C)?\n";

chomp(my $aspect = uc(<STDIN>));

print "Enter the fully qualified name of your associations file:\n";

chomp(my $annotationFile = <STDIN>); 

print "Enter a the fully qualified name of your file with a list of genes for which to find term:\n";

chomp(my $genesFile = <STDIN>);

print "How many genes (roughly) exist within the organism?\n";

chomp(my $totalNum = <STDIN>);

print "Finding terms...\n";

my $ontology   = GO::OntologyProvider::OntologyParser->new(ontologyFile=>$ontologyFile);

my $annotation = GO::AnnotationProvider::AnnotationParser->new(annotationFile=>$annotationFile);

my $termFinder = GO::TermFinder->new(annotationProvider=> $annotation,
				     ontologyProvider  => $ontology,
				     totalNumGenes     => $totalNum,
				     aspect            => $aspect);

my @genes = GenesFromFile($genesFile);

my @pvalues    = $termFinder->findTerms(genes        => \@genes,
					calculateFDR => 1);

# now just print the info back to the client

my $report = GO::TermFinderReport::Text->new();

my $cutoff = 0.05;

my $numHypotheses = $report->print(pvalues  => \@pvalues,
				   numGenes => scalar(@genes),
				   totalNum => $totalNum,
				   cutoff   => $cutoff);

# if they had no significant P-values

if ($numHypotheses == 0){
    
    print "No terms were found for this aspect with a corrected P-value <= $cutoff.\n";
    
}

=pod

=head1 NAME

termFinderClient.pl - interactive client to find significant GO terms for a list of genes

=head1 SYNOPSIS

This program is a very simply client for the GO::TermFinder object,
that prompts a user for the various pieces of information that are
required to determine significant GO terms associated with the list of
genes.  It uses a p-value cut-off (for the corrected p-value) of .05,
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
    -- 1 of 30 --
    GOID    GO:0006790
    TERM    sulfur metabolism
    CORRECTED P-VALUE       3.16574033239377e-23
    UNCORRECTED P-VALUE     4.72498557073697e-25
    NUM_ANNOTATIONS 13 of 19 in the list, vs 51 of 7300 in the genome
    The genes annotated to this node are:
    STR3, MHT1, MET17, MET6, ECM17, MET1, MET2, MET3, MET16, MET28, MET11, MET14, MET10
    
    -- 2 of 30 --
    GOID    GO:0000096
    TERM    sulfur amino acid metabolism
    CORRECTED P-VALUE       2.20564875638675e-21
    UNCORRECTED P-VALUE     3.29201306923395e-23
    NUM_ANNOTATIONS 11 of 19 in the list, vs 29 of 7300 in the genome
    The genes annotated to this node are:
    STR3, MHT1, MET17, MET6, MET1, MET2, MET3, MET16, MET28, MET11, MET14
    
    -- 3 of 30 --
    GOID    GO:0006555
    TERM    methionine metabolism
    CORRECTED P-VALUE       6.3525062302867e-18
    UNCORRECTED P-VALUE     9.48135258251747e-20
    NUM_ANNOTATIONS 9 of 19 in the list, vs 20 of 7300 in the genome
    The genes annotated to this node are:
    STR3, MET17, MET6, MET1, MET2, MET3, MET16, MET11, MET14
    
    -- 4 of 30 --
    GOID    GO:0006520
    TERM    amino acid metabolism
    CORRECTED P-VALUE       8.86722504694052e-16
    UNCORRECTED P-VALUE     1.3234664249165e-17
    NUM_ANNOTATIONS 13 of 19 in the list, vs 175 of 7300 in the genome
    The genes annotated to this node are:
    STR3, MHT1, MET17, MET6, ECM17, MET1, MET2, MET3, MET16, MET28, MET11, MET14, YIL074C
    
    -- 5 of 30 --
    GOID    GO:0006519
    TERM    amino acid and derivative metabolism
    CORRECTED P-VALUE       2.14382262130395e-15
    UNCORRECTED P-VALUE     3.19973525567754e-17
    NUM_ANNOTATIONS 13 of 19 in the list, vs 187 of 7300 in the genome
    The genes annotated to this node are:
    STR3, MHT1, MET17, MET6, ECM17, MET1, MET2, MET3, MET16, MET28, MET11, MET14, YIL074C
    
    -- 6 of 30 --
    GOID    GO:0009308
    TERM    amine metabolism
    CORRECTED P-VALUE       4.57805242229355e-15
    UNCORRECTED P-VALUE     6.8329140631247e-17
    NUM_ANNOTATIONS 13 of 19 in the list, vs 198 of 7300 in the genome
    The genes annotated to this node are:
    STR3, MHT1, MET17, MET6, ECM17, MET1, MET2, MET3, MET16, MET28, MET11, MET14, YIL074C
    
    -- 7 of 30 --
    GOID    GO:0009066
    TERM    aspartate family amino acid metabolism
    CORRECTED P-VALUE       1.64111443018549e-14
    UNCORRECTED P-VALUE     2.44942452266491e-16
    NUM_ANNOTATIONS 9 of 19 in the list, vs 42 of 7300 in the genome
    The genes annotated to this node are:
    STR3, MET17, MET6, MET1, MET2, MET3, MET16, MET11, MET14
    
    -- 8 of 30 --
    GOID    GO:0006082
    TERM    organic acid metabolism
    CORRECTED P-VALUE       1.4994989396095e-13
    UNCORRECTED P-VALUE     2.23805811882014e-15
    NUM_ANNOTATIONS 13 of 19 in the list, vs 258 of 7300 in the genome
    The genes annotated to this node are:
    STR3, MHT1, MET17, MET6, ECM17, MET1, MET2, MET3, MET16, MET28, MET11, MET14, YIL074C
    
    -- 9 of 30 --
    GOID    GO:0019752
    TERM    carboxylic acid metabolism
    CORRECTED P-VALUE       1.4994989396095e-13
    UNCORRECTED P-VALUE     2.23805811882014e-15
    NUM_ANNOTATIONS 13 of 19 in the list, vs 258 of 7300 in the genome
    The genes annotated to this node are:
    STR3, MHT1, MET17, MET6, ECM17, MET1, MET2, MET3, MET16, MET28, MET11, MET14, YIL074C
    
    -- 10 of 30 --
    GOID    GO:0000103
    TERM    sulfate assimilation
    CORRECTED P-VALUE       2.41921570979174e-13
    UNCORRECTED P-VALUE     3.61076971610707e-15
    NUM_ANNOTATIONS 6 of 19 in the list, vs 8 of 7300 in the genome
    The genes annotated to this node are:
    MET14, ECM17, MET1, MET3, MET16, MET10
    
    -- 11 of 30 --
    GOID    GO:0006791
    TERM    sulfur utilization
    CORRECTED P-VALUE       2.41921570979174e-13
    UNCORRECTED P-VALUE     3.61076971610707e-15
    NUM_ANNOTATIONS 6 of 19 in the list, vs 8 of 7300 in the genome
    The genes annotated to this node are:
    MET14, ECM17, MET1, MET3, MET16, MET10
    
    -- 12 of 30 --
    GOID    GO:0000097
    TERM    sulfur amino acid biosynthesis
    CORRECTED P-VALUE       7.24655785462269e-13
    UNCORRECTED P-VALUE     1.08157579919742e-14
    NUM_ANNOTATIONS 6 of 19 in the list, vs 9 of 7300 in the genome
    The genes annotated to this node are:
    STR3, MET17, MET6, MET2, MET28, MET11
    
    -- 13 of 30 --
    GOID    GO:0008652
    TERM    amino acid biosynthesis
    CORRECTED P-VALUE       4.9099398159285e-09
    UNCORRECTED P-VALUE     7.32826838198284e-11
    NUM_ANNOTATIONS 8 of 19 in the list, vs 102 of 7300 in the genome
    The genes annotated to this node are:
    STR3, MET17, MET6, ECM17, MET2, MET28, MET11, YIL074C
    
    -- 14 of 30 --
    GOID    GO:0009309
    TERM    amine biosynthesis
    CORRECTED P-VALUE       8.42568196571993e-09
    UNCORRECTED P-VALUE     1.25756447249551e-10
    NUM_ANNOTATIONS 8 of 19 in the list, vs 109 of 7300 in the genome
    The genes annotated to this node are:
    STR3, MET17, MET6, ECM17, MET2, MET28, MET11, YIL074C
    
    -- 15 of 30 --
    GOID    GO:0000101
    TERM    sulfur amino acid transport
    CORRECTED P-VALUE       9.98458841804645e-06
    UNCORRECTED P-VALUE     1.49023707732037e-07
    NUM_ANNOTATIONS 3 of 19 in the list, vs 5 of 7300 in the genome
    The genes annotated to this node are:
    MMP1, MUP1, SAM3
    
    -- 16 of 30 --
    GOID    GO:0009086
    TERM    methionine biosynthesis
    CORRECTED P-VALUE       9.98458841804645e-06
    UNCORRECTED P-VALUE     1.49023707732037e-07
    NUM_ANNOTATIONS 3 of 19 in the list, vs 5 of 7300 in the genome
    The genes annotated to this node are:
    STR3, MET6, MET2
    
    -- 17 of 30 --
    GOID    GO:0009067
    TERM    aspartate family amino acid biosynthesis
    CORRECTED P-VALUE       0.000797495232457435
    UNCORRECTED P-VALUE     1.19029139172751e-05
    NUM_ANNOTATIONS 3 of 19 in the list, vs 18 of 7300 in the genome
    The genes annotated to this node are:
    STR3, MET6, MET2
    
    -- 18 of 30 --
    GOID    GO:0009069
    TERM    serine family amino acid metabolism
    CORRECTED P-VALUE       0.00129343686399999
    UNCORRECTED P-VALUE     1.93050278208954e-05
    NUM_ANNOTATIONS 3 of 19 in the list, vs 21 of 7300 in the genome
    The genes annotated to this node are:
    MET17, MET2, YIL074C
    
    -- 19 of 30 --
    GOID    GO:0006865
    TERM    amino acid transport
    CORRECTED P-VALUE       0.00520224372539958
    UNCORRECTED P-VALUE     7.76454287373071e-05
    NUM_ANNOTATIONS 3 of 19 in the list, vs 33 of 7300 in the genome
    The genes annotated to this node are:
    MMP1, MUP1, SAM3
    
    -- 20 of 30 --
    GOID    GO:0015837
    TERM    amine/polyamine transport
    CORRECTED P-VALUE       0.0107851312833169
    UNCORRECTED P-VALUE     0.000160972108706222
    NUM_ANNOTATIONS 3 of 19 in the list, vs 42 of 7300 in the genome
    The genes annotated to this node are:
    MMP1, MUP1, SAM3
    
    -- 21 of 30 --
    GOID    GO:0046942
    TERM    carboxylic acid transport
    CORRECTED P-VALUE       0.0141676023767647
    UNCORRECTED P-VALUE     0.00021145675189201
    NUM_ANNOTATIONS 3 of 19 in the list, vs 46 of 7300 in the genome
    The genes annotated to this node are:
    MMP1, MUP1, SAM3
    
    -- 22 of 30 --
    GOID    GO:0015849
    TERM    organic acid transport
    CORRECTED P-VALUE       0.0151086898836281
    UNCORRECTED P-VALUE     0.000225502834084001
    NUM_ANNOTATIONS 3 of 19 in the list, vs 47 of 7300 in the genome
    The genes annotated to this node are:
    MMP1, MUP1, SAM3
    
    -- 23 of 30 --
    GOID    GO:0009070
    TERM    serine family amino acid biosynthesis
    CORRECTED P-VALUE       0.0279454982968357
    UNCORRECTED P-VALUE     0.000417096989505011
    NUM_ANNOTATIONS 2 of 19 in the list, vs 12 of 7300 in the genome
    The genes annotated to this node are:
    MET17, YIL074C

=head1 AUTHORS

Gavin Sherlock, sherlock@genome.stanford.edu

=cut
