use Test;
BEGIN { plan tests => 17 };

# File       : GO-AnnotationProvider-AnnotationParser.t
# Author     : Gavin Sherlock
# Date Begun : March 9th 2002

# $Id$

# This file forms a set of tests for the
# GO::AnnotationProvider::AnnotationParser class

use strict;
use warnings;
use diagnostics;

use GO::OntologyProvider::OntologyParser;

# 'make test' command will be run from one directory up

my $ontologyFile = "t/process.ontology";

my $ontology = GO::OntologyProvider::OntologyParser->new(ontologyFile => $ontologyFile);

# check we're the right type of object

ok($ontology->isa("GO::OntologyProvider::OntologyParser"));

ok($ontology->isa("GO::OntologyProvider"));

# check the object returns a code reference when asked if it can do a
# method that should exist

my @methods = qw (printOntology allNodes rootNode nodeFromId
		  serializeToDisk);

foreach my $method (@methods){

    ok(ref($ontology->can($method)), "CODE");

}

# now we want to check some of the attributes

ok(scalar($ontology->allNodes), 6957); # has been manually checked

my $rootNode = $ontology->rootNode;

# some specifics about the root node

ok($rootNode->goid, "GO:0003673");

ok($rootNode->term, "Gene_Ontology");

ok(scalar($rootNode->childNodes), 1); # should only have 1 child

ok(scalar($rootNode->parentNodes), 0); # should be no parents

# check it's only child is biological_process

ok(($rootNode->childNodes)[0]->goid, "GO:0008150");

ok(($rootNode->childNodes)[0]->term, "biological_process");

# check a random node with plenty of parents

my $node = $ontology->nodeFromId("GO:0042217");

ok(scalar($node->parentNodes), 6);

# now check a random node with lots of children

my $otherNode = $ontology->nodeFromId("GO:0006520");

ok(scalar($otherNode->childNodes), 14);

# now check that each node is valid

my $validNodes = 0;

foreach my $node ($ontology->allNodes){

    $validNodes += $node->isValid;

}

ok(scalar($ontology->allNodes), $validNodes);
