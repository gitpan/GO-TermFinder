use Test;
BEGIN { plan tests => 47 };

# File       : GO-AnnotationProvider-AnnotationParser.t
# Author     : Gavin Sherlock
# Date Begun : March 9th 2002

# $Id$

# This file forms a set of tests for the
# GO::AnnotationProvider::AnnotationParser class

use strict;
use warnings;
use diagnostics;

use GO::AnnotationProvider::AnnotationParser;

# 'make test' command will be run from one directory up

my $associationsFile = "t/gene_association.sgd";

my @methods = qw (nameIsAmbiguous databaseIdsForAmbiguousName
		  ambiguousNames goIdsByDatabaseId goIdsByStandardName
		  goIdsByName standardNameByDatabaseId
		  databaseIdByStandardName databaseIdByName
		  standardNameByName nameIsStandardName
		  nameIsDatabaseId nameIsAnnotated databaseName
		  numAnnotatedGenes allDatabaseIds allStandardNames
		  file serializeToDisk);

my $annotations = GO::AnnotationProvider::AnnotationParser->new(annotationFile => $associationsFile);

# check we're the right kind of object

ok($annotations->isa("GO::AnnotationProvider::AnnotationParser"));

ok($annotations->isa("GO::AnnotationProvider"));

# check the object returns a code reference when asked it it can do a
# method that should exist

foreach my $method (@methods){

    ok(ref($annotations->can($method)), "CODE");

}

# now we want to check some of the attributes

ok($annotations->file, $associationsFile);

ok(scalar($annotations->allDatabaseIds), 6905); # this number has been checked

ok(scalar($annotations->allDatabaseIds), scalar($annotations->allStandardNames));

ok(scalar($annotations->allDatabaseIds), $annotations->numAnnotatedGenes);

# now lets check actin - should be annotated to 5 components

my %components = ("GO:0005884" => undef,       
		  "GO:0000141" => undef,
		  "GO:0000142" => undef,
		  "GO:0005857" => undef,
		  "GO:0000123" => undef);

# NOTE - should think about changing API of goIdsByDatabaseId to
# return an array, rather than a reference to one....

my $goidsRef = $annotations->goIdsByDatabaseId(databaseId => "S0001855",
					       aspect     => 'C');

ok(scalar(@{$goidsRef}), 5);

foreach my $goid (@{$goidsRef}){

    ok(exists($components{$goid})); # should be in the ones we know

}

# now look at TUB1 functions - should be annotated to a single node

my %functions = ("GO:0005200" => undef);

$goidsRef = $annotations->goIdsByDatabaseId(databaseId => "S0004550",
					    aspect     => 'F');

ok(scalar(@{$goidsRef}), 1);

ok(exists($functions{$goidsRef->[0]}));

# now check processes for CDC8 - should be annotated to 6 nodes

my %processes = ("GO:0006233" => undef,
		 "GO:0006235" => undef,
		 "GO:0006261" => undef,
		 "GO:0006276" => undef,
		 "GO:0006280" => undef,
		 "GO:0006281" => undef);

$goidsRef = $annotations->goIdsByDatabaseId(databaseId => "S0003818",
					    aspect     => 'P');

ok(scalar(@{$goidsRef}), 6);

foreach my $goid (@{$goidsRef}){

    ok(exists($processes{$goid})); # should be in the ones we know

}

# now check that all databaseIds are databaseIds

my @databaseIds = $annotations->allDatabaseIds;

my $isDatabaseId = 0;

foreach my $databaseId (@databaseIds){

    $isDatabaseId += $annotations->nameIsDatabaseId($databaseId);

}

ok($isDatabaseId, scalar(@databaseIds));

# now do the same for standard names

my @standardNames = $annotations->allStandardNames;

my $isStandardName = 0;

foreach my $standardName (@standardNames){

    $isStandardName += $annotations->nameIsStandardName($standardName);

}

ok($isStandardName, scalar(@standardNames));

# now check the ambiguous names

my @ambiguousNames = $annotations->ambiguousNames;

my $isAmbiguousName = 0;

foreach my $ambiguousName (@ambiguousNames){

    $isAmbiguousName += $annotations->nameIsAmbiguous($ambiguousName);

}

ok($isAmbiguousName, scalar(@ambiguousNames));

# now check some specific data

ok($annotations->nameIsDatabaseId("S0003818"));

ok($annotations->nameIsStandardName("ACT1"));

ok($annotations->nameIsAmbiguous("HAP1"));

ok(scalar($annotations->databaseIdsForAmbiguousName("HAP1")), 2);

