use Test;
BEGIN { plan tests => 23 };

# File       : GO-Node.t
# Author     : Gavin Sherlock
# Date Begun : March 9th 2002

# $Id$

# This file forms a set of tests for the GO::Node class

use strict;
use warnings;
use diagnostics;

use GO::Node;

my $goid = "GO:0008150";
my $term = "biological_process";

my @methods = qw(addChildNodes addParentNodes addPathToRoot goid term
		 childNodes parentNodes pathsToRoot pathsToAncestor
		 ancestors lengthOfLongestPathToRoot isValid
		 isAParentOf isAChildOf isAnAncestorOf isADescendantOf
		 isLeaf isRoot);

my $node = GO::Node->new(goid => $goid,
			 term => $term);

# check that we're the right kind of object

ok($node->isa("GO::Node"));

# check the object returns a code reference when asked it it can do a
# method that should exist

foreach my $method (@methods){

    ok(ref($node->can($method)), "CODE");

}

# now check attribute values

ok($node->goid, $goid);

ok($node->term, $term);

# now check we get an appropriate error thrown if we miss out a
# required argument

# leave out term

eval {

    $node = GO::Node->new(goid => $goid);

};

ok($@ =~/did not provide a value for the 'term' argument/);

# leave out goid

eval {

    $node = GO::Node->new(term => $term);

};

ok($@ =~/did not provide a value for the 'goid' argument/);
