#!/usr/bin/perl

# $Id: batchGOView.pl,v 1.3 2003/12/17 02:15:36 sherlock Exp $

# Date   : 4th December 2003
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
use diagnostics;
use warnings;

use CGI qw/:all :html3/;

=head1 NAME

batchGOView.pl - batch processor for creating visual output from GO::TermFinder

=head1 SYNOPSIS

batchGoView.pl will read through a number of files, each containing a
list of genes, and will create for each one an html page with a
GO::View, such that you can graphically browse the results.  You need
to provide a .conf file, and then a list of files for processing, each
of which contain a list of genes.  An example .conf file exists in
this directory - edit as appropriate.

Usage:

batchGOView.pl <.conf_file> <file1> <file2> <file3> ... <fileN>

The following usage should give you some output, using the example
files:

batchGOView.pl GoView.conf genes.txt genes2.txt

An html file, batchGOView.html will be created, that will allow you to
browse the results from all of the input lists of genes in a sinple
format.  A frame on the left will have a list of the files that were
input, and the frame on the right will display the results for the
clicked on link.

=cut

use GO::TermFinder;
use GO::AnnotationProvider::AnnotationParser;
use GO::OntologyProvider::OntologyParser;
use GO::View;
use GO::TermFinderReport::Html;
use GO::Utils::File    qw (GenesFromFile);
use GO::Utils::General qw (CategorizeGenes);

$|=1;

##################################################################################
sub Usage{
###################################################################################


    print <<USAGE;

This program takes a list of files, each of which contain a list of
genes, with one gene per line.  It will findTerms for the lists of
genes in each of the GO aspects, and then generate an html page
with a GO::View graphic that summarize the result.

It will use the first supplied argument as the annotation file, the
second argument ontology file, the third argument as the aspect of the
ontology (P, C, or F), and the fourth argument as the expected number
of genes within the organism, and all subsequent files as ones
containing lists of genes.

Usage:

batchGOView.pl <.conf_file> <file1> <file2> <file3> ... <fileN>

e.g.

batchGOView.pl GoView.conf genes.txt genes2.txt

USAGE

    exit;

}

# we need at least 2 arguments, a .conf file, and a file of input
# genes to test

&Usage if (@ARGV < 2);

my $confFile = shift;

my ($annotationFile, $ontologyFile, $aspect, $totalNumGenes, $outDir,
    $geneUrl, $goidUrl, $pvalueCutOff) = &ReadConfFile($confFile);

# now set up the objects we need

my $ontology = GO::OntologyProvider::OntologyParser->new(ontologyFile => $ontologyFile);

my $annotation = GO::AnnotationProvider::AnnotationParser->new(annotationFile=>$annotationFile);

$totalNumGenes ||= $annotation->numAnnotatedGenes;

my $termFinder = GO::TermFinder->new(annotationProvider=> $annotation,
				     ontologyProvider  => $ontology,
				     totalNumGenes     => $totalNumGenes,
				     aspect            => $aspect);

my $report  = GO::TermFinderReport::Html->new();

&GenerateFrameset;

# now open an html file that will have a list of links for all the results

open (LIST, ">".$outDir."batchGOViewList.html")|| die "Cannot create ".$outDir."list.html : $!";

# go through each file

foreach my $file (@ARGV){
   
    print "Analyzing $file\n";

    # get the genes in the file

    my @genes = GenesFromFile($file);

    # now we have to decide which we can analyze
    
    my (@list, @notFound, @ambiguous);

    CategorizeGenes(annotation  => $annotation,
		    genes       => \@genes,
		    ambiguous   => \@ambiguous,
		    unambiguous => \@list,
		    notFound    => \@notFound);

    # technically, those genes that are not found should still
    # contribute to the finding of terms (they will just be considered
    # as unannotated).

    # obviously if no genes are found as unambiguous, then we will
    # find nothing of interest.

    if (!@list){

	print "No known genes exist in $file, so skipping.\n";
	next;

    }

    # now we want to find terms

    my @pvalues = $termFinder->findTerms(genes=>[(@list, @notFound)]);

    # now we hand these off to the GO::View module, to create the image etc.

    my $goView = GO::View->new(-ontologyProvider   => $ontology,
			       -annotationProvider => $annotation,
			       -termFinder         => \@pvalues,
			       -aspect             => $aspect,
			       -configFile         => $confFile,
			       -imageDir           => $outDir,
			       -imageLabel         => "Batch GO::View",
			       -nodeUrl            => $goidUrl,
			       -geneUrl            => $geneUrl,
			       -pvalueCutOff       => $pvalueCutOff);

    # We now want to get the image and map that has hopefully been
    # created by the GO::View module, so we can print it to our html
    # page

    my $imageFile;

    if ($goView->graph) {
	
	$imageFile = $goView->showGraph;
	
    }
    
    my $htmlFile = &GenerateHTMLFile($file, $goView->imageMap, \@pvalues,
				     scalar(@list) + scalar(@notFound), "Terms for $file"); 

    print LIST a({-href=>$htmlFile,
		  -target=>'result'}, $htmlFile), br;

}

close LIST;

sub GenerateHTMLFile{

    my ($file, $map, $pvaluesRef, $numGenes, $title) = @_;

    # work out name of html file
    
    my $htmlFile = $file;

    # delete anything up to and including the last slash

    $htmlFile =~ s/.*\///;

    # delete anything following the last period

    $htmlFile =~ s/\..*//;

    # now add an html suffix

    $htmlFile .= ".html";

    my $fullHtmlFile = $outDir.$htmlFile;

    open (HTML, ">$fullHtmlFile") || die "Cannot create $fullHtmlFile : $!";

    print HTML start_html(-title=>$title);

    print HTML center(h2($title)), hr;

    print HTML $map if defined $map;

    my $numRows = $report->print(pvalues      => $pvaluesRef,
				 aspect       => $aspect,
				 numGenes     => $numGenes,
				 totalNum     => $totalNumGenes,
				 fh           => \*HTML,
				 pvalueCutOff => $pvalueCutOff,
				 geneUrl      => $geneUrl,
				 goidUrl      => $goidUrl);

    if ($numRows == 0){

	print HTML h4(font({-color=>'red'}),
		      center("There were no GO nodes exceeding the p-value cutoff of $pvalueCutOff for the genes in $file."));

    }

    print HTML end_html;

    close HTML;

    return ($htmlFile);

}

sub ReadConfFile{

    my $confFile = shift;

    my %conf;

    open (CONF, $confFile) || die "cannot open $confFile : $!";

    while (<CONF>){

	next if /^\#/; # skip comment lines

	chomp;

	next if /^\s*$/; # skip blank lines, or those without content

	next unless /(.+) = (.+)/;

	my ($param, $value) = ($1, $2);

	$value =~ s/\s$//;

	$conf{$param} = $value;

    }

    close CONF;

    if (!exists $conf{'annotationFile'} || !defined $conf{'annotationFile'}){

	die "Your conf file must specify an annotation file entry.";

    }elsif (!exists $conf{'ontologyFile'} || !defined $conf{'ontologyFile'}){

	die "Your conf file must specify an ontology file entry.";

    }elsif (!exists $conf{'aspect'} || !defined $conf{'aspect'}){

	die "Your conf file must specify an aspect entry.";

    }

    if (!exists $conf{'totalNumGenes'} || !defined $conf{'totalNumGenes'}){

	$conf{'totalNumGenes'} = ""; # simply make it the empty string for now

    }

    if (!exists $conf{'outDir'} || !defined $conf{'outDir'}){

	$conf{'outDir'} = ""; # set to empty string for now

    }

    $conf{'geneUrl'} ||= "";
    $conf{'goidUrl'} ||= "";

    $conf{'pvalueCutOff'} ||= 1;

    # now make sure that file paths are treated relative to the conf file

    my $confDir = "./"; # default

    if ($confFile =~ /(.+)\//){

	$confDir = $1."/"; # adjust if necessary

    }

    foreach my $file ($conf{'annotationFile'}, $conf{'ontologyFile'}, $conf{'outDir'}){

	# $file is an alias for the hash entry

	if ($file !~ /^\//){ # if it's not an absolute path

	    $file = $confDir.$file; # add the confDir on the front

	}

    }

    # return a hash slice

    return @conf{qw(annotationFile ontologyFile aspect totalNumGenes
		    outDir geneUrl goidUrl pvalueCutOff)};

}

sub GenerateFrameset{

# start an index file that a user can use to browse the output data,
# using frames

    open (FRAMES, ">".$outDir."batchGOView.html") || die "Cannot create ${outDir}batchGOView.html : $!";

    print FRAMES frameset({-cols         => "100, *",
			   -marginheight => '0',
			   -marginwidth  => '0',
			   -frameborder  => '1',
			   -border       => '1'},
			  
			  frame({'-name'       => "list",
				 -src          => "batchGOViewList.html",
				 -marginwidth  => 0,
				 -marginheight => 0,
				 -border       => 1}),
		   
			  frame({'-name'       =>'result',
				 -marginwidth  => 0,
				 -marginheight => 0,
				 -border       => 1}));

    close FRAMES;

}

=pod

=head1 AUTHOR

Gavin Sherlock

sherlock@genome.stanford.edu

=cut
