#!/usr/bin/perl

# $Id: analyze.pl,v 1.3 2003/10/17 23:28:33 sherlock Exp $

# Date   : 16th October 2003
# Author : Gavin Sherlock

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

$|=1;

###################################################################################
sub Usage{
###################################################################################


    print <<USAGE;

This program takes a list of files, each of which contain a list of
genes, with one gene per line.  It will findTerms for the lists of
genes in each of the GO aspects, outputting the results to a file
named for the original file, but with a .terms extension.  It will only
output terms with a corrected P-value of <= 0.05.

It will use the first supplied argument as the annotation file, the
second argument as the expected number of genes within the organism,
and all subsequent files as ones containing lists of genes.  You need
to provide the ontology files in the same directory from which you are
executing this script, and the GO-TermFinder libraries must be in your
path.

Usage:

analyze.pl <annotation_file> <numGenes> <file1> <file2> <file3> ... <fileN>

e.g.

analyze.pl gene_association.sgd 7200 file1.txt file2.txt file3.txt

USAGE

    exit;

}

# we need at least 3 arguments, an annotation file, the number of
# genes in the genome, and a file of input genes to test

&Usage if (@ARGV < 3);

# now get our annotation file and number of genes

my $annotationFile = shift;
my $numGenes       = shift;

# now set up the objects we need

my $process   = GO::OntologyProvider::OntologyParser->new(ontologyFile => "process.ontology");
my $component = GO::OntologyProvider::OntologyParser->new(ontologyFile => "component.ontology");
my $function  = GO::OntologyProvider::OntologyParser->new(ontologyFile => "function.ontology");

my $annotation = GO::AnnotationProvider::AnnotationParser->new(annotationFile=>$annotationFile);

my $termFinderP = GO::TermFinder->new(annotationProvider=> $annotation,
				      ontologyProvider  => $process,
				      totalNumGenes     => $numGenes,
				      aspect            => 'P');


my $termFinderC = GO::TermFinder->new(annotationProvider=> $annotation,
				      ontologyProvider  => $component,
				      totalNumGenes     => $numGenes,
				      aspect            => 'C');

my $termFinderF = GO::TermFinder->new(annotationProvider=> $annotation,
				      ontologyProvider  => $function,
				      totalNumGenes     => $numGenes,
				      aspect            => 'F');


# now go through each file

foreach my $file (@ARGV){

    print "Analyzing $file\n";

    my $outfile = $file.".terms";

    open (IN,  $file)        || die "Cannot open $file : $!";
    open (OUT, ">".$outfile) || die "Cannot make $outfile : $!"; 

    my @genes;
    my @lines;

    my $var = chr(13); # to deal with Mac end of line 

    while (<IN>){

	if (/$var/o){ 

	    # if it's a Mac file multiple lines get read at once, so
	    # we have to split on the end-of line character

	    @lines = split($var, $_);

	}else{

	    @lines = ($_);

	}

	foreach my $gene (@lines){

	    $gene =~ s/\cM//g; # remove Control-M characters
	    
	    $gene =~ s/\s+$//; # remove any trailing or leading whitespace
	    $gene =~ s/^\s//;
	    
	    next unless $gene;
	    
	    push (@genes, $gene);

	}

    }

    close IN;

    my (@list, @notFound, @ambiguous);

    foreach my $gene (@genes) {

	if ($annotation->nameIsAmbiguous($gene)){

	    push(@ambiguous, $gene);
	    print "$gene is ambiguous\n";
	    next;

	}

	my $name = $annotation->standardNameByName($gene);

	if (defined $name){

	    push(@list, $gene);

	}else{

	    push(@notFound, $gene);

	}

    }

    if (@list){

	print OUT "The following gene(s) will be considered:\n\n";

	foreach my $gene (@list){

	    print OUT $gene, "\t", $annotation->standardNameByName($gene), "\n";

	}

	print OUT "\n";

    }else{

	print OUT "None of the gene names were recognized\n";
	print OUT "They were:\n\n";
	print OUT join("\n", @notFound), "\n";
	close OUT;

	next;

    }

    if (@ambiguous){

	print OUT "The following gene(s) are ambiguously named, and so will not be used:\n";
	print OUT join("\n", @ambiguous), "\n\n";

    }

    if (@notFound){

	print OUT "The following gene(s) were not recognized, and will not be considered:\n\n";
	print OUT join("\n", @notFound), "\n\n";

    }

    foreach my $termFinder ($termFinderP, $termFinderC, $termFinderF){

	print OUT "Finding terms for ", $termFinder->aspect, "\n\n";

	my @pvalues = $termFinder->findTerms(genes=>\@list);

	my $hypothesis = 1;
	
	foreach my $pvalue (@pvalues){
	    
	    next if ($pvalue->{CORRECTED_PVALUE} > 0.05);
	
	    print OUT "-- $hypothesis of ", scalar @pvalues, " --\n",
	    
	    "GOID\t", $pvalue->{NODE}->goid, "\n",
	    
	    "TERM\t", $pvalue->{NODE}->term, "\n",
	    
	    "CORRECTED P-VALUE\t", $pvalue->{CORRECTED_PVALUE}, "\n",
	    
	    "UNCORRECTED P-VALUE\t", $pvalue->{PVALUE}, "\n",

	    "NUM_ANNOTATIONS\t", $pvalue->{NUM_ANNOTATIONS}, " of ", scalar (@list), " in the list, vs ", $pvalue->{TOTAL_NUM_ANNOTATIONS}, " of $numGenes in the genome\n",

	    "The genes annotated to this node are:\n",

	    join(", ", values(%{$pvalue->{ANNOTATED_GENES}})), "\n";
	    
	    $hypothesis++;
	    
	}

	# if they had no significant P-values

	if ($hypothesis == 1){

	    print OUT "No terms were found for this aspect with a corrected P-value <= 0.05.\n";

	}

	print OUT "\n\n";

    }

    close OUT;    
    
}

