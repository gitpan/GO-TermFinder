#!/usr/bin/perl

# $Id: children.pl,v 1.4 2003/11/26 18:46:47 sherlock Exp $

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

use GO::OntologyProvider::OntologyParser;

my ($goid, $ontologyFile) = @ARGV;

&Usage("You must provide a goid")           if (!$goid);
&Usage("You must provide an ontology file") if (!$ontologyFile);
&Usage("Your ontology file does not exist") if (!-e $ontologyFile);

my $ontology = GO::OntologyProvider::OntologyParser->new(ontologyFile => $ontologyFile);

my $node = $ontology->nodeFromId($goid);

my @children = $node->childNodes;

print "\n";

if (@children){

    print "Children of $goid (", $node->term, ") : \n\n";

    foreach my $child (@children){

	print $child->goid, "    ", $child->term, "\n";

    }

}else{

    print $goid, " has no children.\n";

}

print "\n";

sub Usage{

    my $message = shift;

    print $message, ".\n\n";

    print "Usage :

ancestors.pl <goid> <ontology_file>\n\n";

    exit;

}

=pod

=head1 NAME

children.pl - prints children of a supplied GO node

=head1 SYNOPSIS

children.pl simply takes as input a GOID, and an ontology file, and prints out
the children of that GO node, e.g.:

    >children.pl GO:0008150 process.ontology

    Children of GO:0008150 (biological_process) : 

    GO:0000004    biological_process unknown
    GO:0016032    viral life cycle
    GO:0007582    physiological processes
    GO:0007275    development
    GO:0009987    cellular process
    GO:0008371    obsolete
    GO:0007610    behavior

=head1 AUTHORS

Gavin Sherlock, sherlock@genome.stanford.edu

=cut
