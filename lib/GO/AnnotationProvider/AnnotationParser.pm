package GO::AnnotationProvider::AnnotationParser;

# File       : AnnotationParser.pm
# Authors    : Elizabeth Boyle; Gavin Sherlock
# Date Begun : Summer 2001
# Rewritten  : September 25th 2002

# $Id: AnnotationParser.pm,v 1.22 2003/03/03 18:11:45 sherlock Exp $

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

=pod

=head1 NAME

GO:AnnotationProvider::AnnotationParser

=head1 SYNOPSIS

GO:AnnotationProvider::AnnotationParser - reads an Gene Ontology
gene associations file, and provides methods by which to retrieve the
GO annotations for the an annotated entity.

    my $annotationParser = GO::AnnotationProvider::AnnotationParser->new(annotationFile => "data/gene_association.sgd");

    my $geneName = "AAT2";

    print "GO associations for gene: ", join (" ", $annotationParser->goIdsByName(name   => $geneName,
										  aspect => 'P')), "\n";

    print "Database ID for gene: ", $annotationParser->databaseIdByName($geneName), "\n";

    print "Database name: ", $annotationParser->databaseName(), "\n";

    print "StandardName gene name: ", $annotationParser->standardNameForName($geneName), "\n";

    my $i;

    my @geneNames = $annotationParser->allStandardNames();

    foreach $i (0..10) {

        print "$geneNames[$i]\n";

    }

=head1 DESCRIPTION

GO::AnnotationProvider::AnnotationParser is a concrete subclass of
GO::AnnotationProvider, and creates a data structure mapping gene
names to GO annotations by parsing a file of annotations provided by
the Gene Ontology Consortium.

This package provides object methods for retrieving GO annotations
that have been parsed from a 'gene associations' file, provided by
the gene ontology consortium.  The format for the file is:

Lines beginning with a '!' character are comment lines.

    Column  Cardinality   Contents          
    ------  -----------   -------------------------------------------------------------
        0       1         Database abbreviation for the source of annotation (eg SGD)
        1       1         Database identifier of the annotated entity
        2       1         Standard name of the annotated entity
        3       0,1       NOT (if a gene is specifically NOT annotated to the term)
        4       1         GOID of the annotation     
        5       1,n       Reference(s) for the annotation 
        6       1         Evidence code for the annotation
        7       0,n       With or From (a bit mysterious)
        8       1         Aspect of the Annotation (C, F, P)
        9       0,1       Name of the product being annotated
       10       0,n       Alias(es) of the annotated product
       11       1         type of annotated entity (one of gene, transcript, protein)
       12       1,2       taxonomic id of the organism encoding and/or using the product
       13       1         Date of annotation YYYYMMDD
       14       1         Assigned_by : The database which made the annotation

Columns are separated by tabs.  For those entries with a cardinality
greater than 1, multiple entries are pipe , |, delimited.

Further details can be found at:

http://www.geneontology.org/doc/GO.annotation.html#file

The following assumptions about the file are made (and should be true):

    1.  All aliases appear for all entries of a given annotated product
    2.  The database identifiers are unique, in that two different
        entities cannot have the same database id.

=head1 TODO

Also see the TODO list in the parent, GO::AnnotationProvider.

 1.  Add in methods that will allow retrieval of evidence codes with
     the annotations for a particular entity.

 2.  Add in methods that return all the annotated entities for a
     particular GOID.

 3.  Add in the ability to request only annotations either including
     or excluding particular evidence codes.  Such evidence codes
     could be provided as an anonymous array as the value of a named
     argument.

 4.  Same as number 3, except allow the retrieval of annotated
     entities for a particular GOID, based on inclusion or exclusion
     of certain evidence codes.

 These first four items will require a reworking of how data are
 stored on the backend, and thus the parsing code itself, though it
 should not affect any of the already existing API.

 5.  Instead of 'use'ing Storable, 'require' it instead, only at the
     point of use, which will mean that AnnotationParser can be
     happily used in the absence of Storable, just without those
     functions that need it.

 6.  Extend the ValidateFile class method to check that an entity
     should never be annotated to the same node twice, with the same
     evidence, with the same reference.

 7.  An additional checker, that uses an AnnotationProvider in
     conjunction with an OntologyProvider, would be useful, that
     checks that some of the annotations themselves are valid, ie
     that no entities are annotated to the 'unknown' node in a
     particular aspect, and also to another node within that same
     aspect.  Can annotations be redundant? ie, if an entity is
     annotated to a node, and an ancestor of the node, is that
     annotation redundant?  Does it depend on the evidence codes and
     references.  Or are such annotations reinforcing?  These things
     are useful to consider when formulating the confidence which can
     be attributed to an annotation.

=cut

use strict;
use warnings;
use diagnostics;

use Storable qw (nstore);

use vars qw (@ISA $PACKAGE $VERSION);

use GO::AnnotationProvider;
@ISA = qw (GO::AnnotationProvider);

$PACKAGE = "GO::AnnotationProvider::AnnotationParser";
$VERSION = "0.1";

# CLASS Attributes
#
# These should be considered as constants, and are initialized here

my $DEBUG = 0;

# constants for instance attribute name

my $kDatabaseName      = $PACKAGE.'::__databaseName';      # stores the name of the annotating database
my $kFileName          = $PACKAGE.'::__fileName';          # stores the name of the file used to instantiate the object
my $kNameToIdMap       = $PACKAGE.'::__nameToIdMap';       # stores a map of all names for a gene to the database id
my $kAmbiguousNames    = $PACKAGE.'::__ambiguousNames';    # stores the database id's for all ambiguous names
my $kIdToStandardName  = $PACKAGE.'::__idToStandardName';  # stores a map of database id's to standard names of all entities
my $kStandardNameToId  = $PACKAGE.'::__StandardNameToId';  # stores a map of standard names to their database id's
my $kGoids             = $PACKAGE.'::__goids';             # stores all the goid annotations
my $kNumAnnotatedGenes = $PACKAGE.'::__numAnnotatedGenes'; # stores number of genes with annotations, per aspect

my $kTotalNumAnnotatedGenes = $PACKAGE.'::__totalNumAnnotatedGenes'; # total number of annotated genes

# constants to describe what is in which column in the annotation file

my $kDatabaseNameColumn = 0;
my $kDatabaseIdColumn   = 1;
my $kStandardNameColumn = 2;
my $kNotColumn          = 3;
my $kGoidColumn         = 4;
my $kReferenceColumn    = 5;
my $kEvidenceColumn     = 6;
my $kWithColumn         = 7;
my $kAspectColumn       = 8;
my $kNameColumn         = 9;
my $kAliasesColumn      = 10;
my $kEntityTypeColumn   = 11;
my $kTaxonomicIDColumn  = 12;
my $kDateColumn         = 13;
my $kAssignedByColumn   = 14;

# the following hash of anonymous arrays indicates for each column
# what the maximum and minimum number of entries per column can be.
# If no maximum is indicated, then the maximum is equal to the
# minimum, and exactly that number of entries must exist.

my %kColumnsToCardinality = ($kDatabaseNameColumn => [1     ],
			     $kDatabaseIdColumn   => [1     ],
			     $kStandardNameColumn => [1     ],
			     $kNotColumn          => [0,   1],
			     $kGoidColumn         => [1     ],
			     $kReferenceColumn    => [1, "n"],
			     $kEvidenceColumn     => [1     ],
			     $kWithColumn         => [0, "n"],
			     $kAspectColumn       => [1     ],
			     $kNameColumn         => [0,   1],
			     $kAliasesColumn      => [0, "n"],
			     $kEntityTypeColumn   => [1     ],
			     $kTaxonomicIDColumn  => [1,   2],
			     $kDateColumn         => [1     ],
			     $kAssignedByColumn   => [1     ]);

my $kNumColumnsInFile = scalar keys %kColumnsToCardinality;

############################################################################
#
# CLASS METHODS
#
############################################################################

# Public

############################################################################
sub Usage{
############################################################################
# This class method simply prints out a usage statement, along with an error 
# message, if one was passed in.
#
# Usage :
#
#    GO::AnnotationProvider::AnnotationParser->Usage(), "\n";
#

    my ($class, $message) = @_;

    defined $message && print $message."\n\n";

    print 'The constructor expects one of two arguments, either a
\'annotationFile\' argument, or and \'objectFile\' argument.  When
instantiated with an annotationFile argument, it expects it to
correspond to an annotation file created by one of the GO consortium
members, according to their file format.  When instantiated with an
objectFile argument, it expects to open a previously created
annotationParser object that has been serialized to disk (see the
serializeToDisk method).

Usage:

    my $annotationParser = '.$PACKAGE.'->new(annotationFile => $file);

    my $annotationParser = '.$PACKAGE.'->new(objectFile => $file);
';

}

############################################################################
sub ValidateFile{
############################################################################
# This class method reads an annotation file, and returns a reference
# to an array of errors that are present within the file.  The errors
# are simply strings, each beginning with "Line $lineNo" where $lineNo
# is the number of the line in the file where the error was found.
#
# Usage:
#
#    my $errorsRef = GO::AnnotationProvider::AnnotationParser->ValidateFile(annotationFile => $file);

    my ($class, %args) = @_;
    
    my $file = $args{'annotationFile'} || $class->_handleMissingArgument(argument => 'annotationFile');
    
    open (ANNOTATIONS, $file) || die "$PACKAGE cannot open $file : $!";
    
    my (@errors, @line);

    my ($databaseId, $standardName, $aliases);
    my (%idToName, %idToAliases);

    my $lineNo = 0;
    
    while (<ANNOTATIONS>){
	
	++$lineNo;
	
	next if /^!/; # skip comment lines
	
	chomp;

	next unless $_; # skip an empty line, if there is one
	
	@line = split("\t", $_, -1);
	
	if (scalar @line != $kNumColumnsInFile){ # doesn't have the correct number of columns
	    
	    push (@errors, "Line $lineNo has ". scalar @line. "columns, instead of $kNumColumnsInFile.");
	    
	}
	
	$class->__CheckCardinalityOfColumns(\@errors, \@line, $lineNo);
	
	# now want to deal with sanity checks...
	
	($databaseId, $standardName, $aliases) = @line[$kDatabaseIdColumn, $kStandardNameColumn, $kAliasesColumn];
	
	next if ($databaseId eq ""); # will have given incorrect cardinality, but nothing more we can do with it

	if (!exists $idToName{$databaseId}){

	    $idToName{$databaseId} = $standardName;

	}elsif ($idToName{$databaseId} ne $standardName){

	    push (@errors, "Line $lineNo : $databaseId has more than one standard name : $idToName{$databaseId} and $standardName.");

	}

	if (!exists $idToAliases{$databaseId}){

	    $idToAliases{$databaseId} = $aliases;

	}elsif($idToAliases{$databaseId} ne $aliases){

	    push (@errors, "Line $lineNo : $databaseId has more than one collections of aliases : $idToAliases{$databaseId} and $aliases.");

	}   

    }
    
    close ANNOTATIONS || die "$PACKAGE cannot close $file : $!";
    
    return \@errors;

}

############################################################################
sub __CheckCardinalityOfColumns{
############################################################################
# This method checks the cardinality of each column on a line
#
# Usage:
#
#    $class->__CheckCardinalityOfColumns(\@errors, \@line, $lineNo);
    
    my ($class, $errorsRef, $lineRef, $lineNo) = @_;
    
    my ($cardinality, $min, $max);
    
    foreach my $column (sort {$a<=>$b} keys %kColumnsToCardinality){
	
	($min, $max) = @{$kColumnsToCardinality{$column}}[0,1];
	
	$cardinality = $class->__GetCardinality($lineRef->[$column], $errorsRef, $lineNo);
	
	if (!defined $max){ # must have a defined number of entries

	    if ($cardinality != $min){
		
		push (@{$errorsRef}, "Line $lineNo : column $column has a cardinality of $cardinality, instead of $min.");
		
	    }
	    
	}else{ # there's a range of allowed number of entries
	    
	    if ($cardinality < $min){ # check if less than minimum
		
		push (@{$errorsRef}, "Line $lineNo : column $column has a cardinality of $cardinality, which is less than the required $min.");
		
	    }elsif ($kColumnsToCardinality{$column}->[1] ne 'n' &&
		    $cardinality > $max){ # check if more than maximum
		
		push (@{$errorsRef}, "Line $lineNo : column $column has a cardinality of $cardinality, which is more than the allowed $max.");
		
	    }
	    
	}
	
    }

}

############################################################################
sub __GetCardinality{
############################################################################
# This private method returns an integer that indicates the
# cardinality of a text string, where multiple entries are assumed to
# be seperated by the pipe character (|).  In addition, it checks
# whether there are null or whitespace only entries.
#
# Usage:
#
#    my $cardinality = $class->__GetCardinality($string);

    my ($class, $string, $errorsRef, $lineNo) = @_;

    my $cardinality;

    if (!defined $string || $string eq ""){

	$cardinality = 0;

    }else{

	my @entries = split(/\|/, $string, -1);

	foreach my $entry (@entries){

	    if (!defined $entry){

		push (@{$errorsRef}, "Line $lineNo : There is an undefined value in the string $string.");

	    }elsif ($entry =~ /^\s+$/){

		push (@{$errorsRef}, "Line $lineNo : There is a white-space only value in the string $string.");

	    }

	}

	$cardinality = scalar @entries;

    }

    return $cardinality;

}

############################################################################
#
# Constructor, and initialization methods.
#
# All initialization methods are private, except, of course, for the
# new() method.
#
############################################################################

############################################################################
sub new{
############################################################################
# This is the constructor for an AnnotationParser object.
#
# The constructor expects one of two arguments, either a
# 'annotationFile' argument, or an 'objectFile' argument.  When
# instantiated with an annotationFile argument, it expects it to
# correspond to an annotation file created by one of the GO consortium
# members, according to their file format.  When instantiated with an
# objectFile argument, it expects to open a previously created
# annotationParser object that has been serialized to disk (see the
# serializeToDisk method).
#
# Usage:
#
#    my $annotationParser = GO::AnnotationProvider::AnnotationParser->new(annotationFile => $file);
#
#    my $annotationParser = GO::AnnotationProvider::AnnotationParser->new(objectFile => $file);

    my ($class, %args) = @_;

    my $self;

    if (exists($args{'annotationFile'})){

	$self = {};

	bless $self, $class;

	$self->__init($args{'annotationFile'});

    }elsif (exists($args{'objectFile'})){

	$self = Storable::retrieve($args{'objectFile'}) || die "Could not instantiate $PACKAGE object from objectFile : $!";

	$self->__setFile($args{'objectFile'});

    }else{

	$class->Usage("An annotationFile or objectFile argument must be provided.");
	die;

    }

    return ($self);

}

############################################################################
sub __init{
############################################################################
# This private method initializes the object by reading in the data
# from the annotation file.
#
# Usage :
#
#    $self->__init($file);
#

    my ($self, $file) = @_;

    $self->__setFile($file);

    open (ANNOTATIONS, $file) || die "$PACKAGE cannot open $file : $!";

    # now read through annotations file

    my (@line, $databaseId, $goid, $aspect, $standardName, $aliases);
    
    while (<ANNOTATIONS>){

	next if /^!/; # skip commented lines

	chomp;

	next unless $_; # skip an empty line, if there is one

	@line = split("\t", $_, -1);

	next if $line[$kNotColumn] eq 'NOT'; # skip annotations NOT to a GOID

	($databaseId, $goid, $aspect) = @line[$kDatabaseIdColumn, $kGoidColumn, $kAspectColumn];
	($standardName, $aliases)     = @line[$kStandardNameColumn, $kAliasesColumn];

	if ($databaseId eq ""){

	    print "On line $. there is a missing databaseId, so it will be ignored.\n";
	    next;

	}

	# record the source of the annotation
	
	$self->{$kDatabaseName} = $line[$kDatabaseNameColumn] if (!exists($self->{$kDatabaseName}));

	# now simply store the GOID
	
	$self->__storeGOID($databaseId, $goid, $aspect);
	    
	# now map the standard name and all aliases to the database id
	    
	$self->__mapNamesToDatabaseId($databaseId, $standardName, $aliases);

    }

    close (ANNOTATIONS) || die "AnnotationParser can't close $file: $!";

    # now count up how many annotated things we have

    foreach my $databaseId (keys %{$self->{$kGoids}}){

	$self->{$kTotalNumAnnotatedGenes}++;

	foreach my $aspect (keys %{$self->{$kGoids}{$databaseId}}){

	    $self->{$kNumAnnotatedGenes}{$aspect}++;

	}

    }

}

############################################################################
sub __setFile{
############################################################################
# This method sets the name of the file used for construction.
# 
# Usage:
#
#    $self->__setFile($file);
#

    my ($self, $file) = @_;

    $self->{$kFileName} = $file;

}

############################################################################
sub __mapNamesToDatabaseId{
############################################################################
# This private method maps all names and aliases to the databaseId of
# an entity.  It also maps the databaseId to itself, to facilitate a
# single way of mapping any identifier to the database id.
#
# Usage :
#
#    $self->__mapNamesToDatabaseId($databaseId, $standardName, $aliases);
#
# where $aliases is a pipe-delimited list of aliases

    my ($self, $databaseId, $standardName, $aliases) = @_;    

    return if (exists $self->{$kIdToStandardName}{$databaseId}); # we've already seen this databaseId

    my @aliases = split(/\|/, $aliases);

    my %seen; # sometimes an alias will be the same as the standard name

    foreach my $name ($databaseId, $standardName, @aliases){

	next if exists ($seen{$name});

	if (exists $self->{$kNameToIdMap}{$name}){ # it must be ambiguous

	    # so record what it maps to

	    push (@{$self->{$kAmbiguousNames}{$name}}, $databaseId);                   # current databaseId
	    push (@{$self->{$kAmbiguousNames}{$name}}, $self->{$kNameToIdMap}{$name}); # previous databaseId

	    # and delete the old one from the unambiguous mapping

	    delete $self->{$kNameToIdMap}{$name};

	}elsif (exists $self->{$kAmbiguousNames}{$name}){ # we already know it's ambiguous

	    push (@{$self->{$kAmbiguousNames}{$name}}, $databaseId);

	}else{ # otherwise simply map it unambiguously for now

	    $self->{$kNameToIdMap}{$name} = $databaseId;

	}

	$seen{$name} = undef; # remember that we've seen the name

    }

    $self->{$kIdToStandardName}{$databaseId}   = $standardName; # record the standard name for the database id
    $self->{$kStandardNameToId}{$standardName} = $databaseId;   # also make the reverse look up

}

############################################################################
sub __storeGOID{
############################################################################
# This private method stores a GOID for an entity, on a per
# aspect basis, in hash.
#
# Usage:
#
#    $self->__storeGOID($databaseId, $goid, $aspect);
#

    my ($self, $databaseId, $goid, $aspect) = @_;

    $self->{$kGoids}{$databaseId}{$aspect}{$goid} = undef;

}

##############################################################################
#
# PUBLIC INSTANCE METHODS
#
##############################################################################

# some methods dealing with ambiguous names

##############################################################################
sub nameIsAmbiguous{
##############################################################################
# This public method returns a boolean to indicate whether a name is
# ambiguous, ie whether the name might map to more than one entity
# (and therefore more than one databaseId).
#
# Usage:
#
#    if ($annotationParser->nameIsAmbiguous($name)){
#
#               # do something useful....or not....
#
#    }

    my ($self, $name) = @_;

    die "You must supply a name to nameIsAmbiguous" if !defined ($name);

    return exists $self->{$kAmbiguousNames}{$name};

}

############################################################################
sub databaseIdsForAmbiguousName{
############################################################################
# This public method returns an array of database identifiers for an
# ambiguous name.  If the name is not ambiguous, or not recognized, an
# empty list will be returned.
#
# Usage:
#
#    my @databaseIds = $annotationParser->databaseIdsForAmbiguousName($name);

    my ($self, $name) = @_;

    die "You must supply a name to databaseIdsForAmbiguousName" if !defined ($name);
    
    if ($self->nameIsAmbiguous($name)){
	
	return @{$self->{$kAmbiguousNames}{$name}};

    }else{

	return ();

    }

}

############################################################################
sub ambiguousNames{
############################################################################
# This method returns an array of names, which from the annotation file have
# been deemed to be ambiguous.
#
# Usage:
#
#    my @ambiguousNames = $annotationParser->ambiguousNames();

    return keys %{$_[0]->{$kAmbiguousNames}};

}

# methods for retrieving GO annotations for entities

############################################################################
sub goIdsByDatabaseId{
############################################################################
# This public method returns a reference to an array of GOIDs that are associated
# with the supplied databaseId for a specific aspect.  If no
# annotations are associated with that databaseId in that aspect, then
# a reference to an empty array will be returned.  If the databaseId is not
# recognized, then undef will be returned.
#
# Usage:
#
#    my $goidsRef = $annotationParser->goIdsByDatabaseId(databaseId => $databaseId,
#                                                        aspect     => <P|F|C>);
#

    my ($self, %args) = @_;

    my $aspect     = $args{'aspect'}     || $self->_handleMissingArgument(argument => 'aspect');
    my $databaseId = $args{'databaseId'} || $self->_handleMissingArgument(argument => 'databaseId');

    if (exists $self->{$kGoids}{$databaseId}){ # it's a recognised database id

	if (exists $self->{$kGoids}{$databaseId}{$aspect}){ # it has annotations

	    return [keys %{$self->{$kGoids}{$databaseId}{$aspect}}];

	}else{ # it has no annotations
	    
	    return []; # reference to empty array

	}

    }else{ # it's not a recognized database id

	return undef;

    }

}

############################################################################
sub goIdsByStandardName{
############################################################################
# This public method returns a reference to an array of GOIDs that are
# associated with the supplied standardName for a specific aspect.  If
# no annotations are associated with the entity with that standard
# name in that aspect, then a a reference to an empty list will be
# returned.  If the supplied name is not used as a standard name, then
# undef will be returned.
#
# Usage:
#
#    my $goidsRef = $annotationParser->goIdsByStandardName(standardName => $databaseId,
#                                                          aspect       => <P|F|C>);
#

    my ($self, %args) = @_;

    my $aspect       = $args{'aspect'}       || $self->_handleMissingArgument(argument => 'aspect');
    my $standardName = $args{'standardName'} || $self->_handleMissingArgument(argument => 'standardName');

    my $databaseId = $self->databaseIdByStandardName($standardName);

    return undef if !defined $databaseId; # there is no such standard name

    return $self->goIdsByDatabaseId(databaseId => $databaseId,
				    aspect     => $aspect);

}

############################################################################
sub goIdsByName{
############################################################################
# This public method returns a reference to an array of GO IDs that
# are associated with the supplied name for a specific aspect.  If
# there are no GO associations for the entity corresponding to the
# supplied name in the provided aspect, then a reference to an empty
# list will be returned.  If the supplied name does not correspond to
# any entity, then undef will be returned.  Because the name can be
# any of the databaseId, the standard name, or any of the aliases, it
# is possible that the name might be ambiguous.  Clients of this
# object should first test whether the name they are using is
# ambiguous, using the nameIsAmbiguous() method, and handle it
# accordingly.  If an ambiguous name is supplied, then it will die.
#
# Usage:
#
# my $goidsRef = $annotationParser->goIdsByName(name   => $name,
#                                               aspect => <P|F|C>);
#

    my ($self, %args) = @_;

    my $aspect = $args{'aspect'} || $self->_handleMissingArgument(argument => 'aspect');
    my $name   = $args{'name'}   || $self->_handleMissingArgument(argument => 'name'); 

    die "You have supplied an ambiguous name to goIdsByName" if ($self->nameIsAmbiguous($name));

    my $databaseId = $self->databaseIdByName($name);

    return undef if !defined $databaseId; # there is no such name

    return $self->goIdsByDatabaseId(databaseId => $databaseId,
				    aspect     => $aspect);

}

# methods for mapping different types of name to each other

############################################################################
sub standardNameByDatabaseId{
############################################################################
# This method returns the standard name for a database id
#
# Usage:
#
#    my $standardName = $annotationParser->standardNameByDatabaseId($databaseId);
#

    my ($self, $databaseId) = @_;

    die "You must supply a databaseId to standardNameByDatabaseId" if !defined ($databaseId);

    return ($self->{$kIdToStandardName}{$databaseId});

}

############################################################################
sub databaseIdByStandardName{
############################################################################
# This method returns the database id for a standard name
#
# Usage:
#
#    my $databaseId = $annotationParser->databaseIdByStandardName($standardName);
#

    my ($self, $standardName) = @_;

    die "You must supply a standardName to databaseIdByStandardName" if !defined ($standardName);

    return ($self->{$kStandardNameToId}{$standardName});

}

############################################################################
sub databaseIdByName{
############################################################################
# This method returns the database id for any identifier for a gene
# (eg by databaseId itself, by standard name, or by alias).  If the
# used name is ambiguous, then the program will die.  Thus clients
# should call the nameIsAmbiguous() method, prior to using this
# method.  If the name does not map to any databaseId, then undef will
# be returned.
#
# Usage:
#
#    my $databaseId = $annotationParser->databaseIdByName($name);
#

    my ($self, $name) = @_;

    die "You must supply a name to databaseIdByName" if !defined ($name);

    die "You have supplied an ambiguous name to databaseIdByName" if ($self->nameIsAmbiguous($name));

    return ($self->{$kNameToIdMap}{$name});

}

############################################################################
sub standardNameByName{
############################################################################
# This public method returns the standard name for the the gene
# specified by the given name.  Because a name may be ambiguous, the
# nameIsAmbiguous() method should be called first.  If an ambiguous
# name is supplied, then it will die with an appropriate error
# message.  If the name does not map to a standard name, then undef
# will be returned.
# 
# Usage:
#
#    my $standardName = $annotationParser->standardNameByName($name);


    my ($self, $name) = @_;

    die "You must supply a name to standardNameByName" if !defined ($name);

    die "You have supplied an ambiguous name to standardNameByName" if ($self->nameIsAmbiguous($name));

    return $self->{$kIdToStandardName}{$self->databaseIdByName($name)};

}

# other methods relating to names

############################################################################
sub nameIsStandardName{
############################################################################
# This method returns a boolean to indicate whether the supplied name is
# used as a standard name.
#
# Usage :
#
# if ($annotationParser->nameIsStandardName($name)){
#
#     # do something
#
# }



    my ($self, $name) = @_;

    die "You must supply a name to nameIsStandardName" if !defined($name);

    return exists ($self->{$kStandardNameToId}{$name});

}

############################################################################
sub nameIsDatabaseId{
############################################################################
# This method returns a boolean to indicate whether the supplied name is
# used as a standard name.
#
# Usage :
#
# if ($annotationParser->nameIsDatabaseId($name)){
#
#     # do something
#
# }

    my ($self, $databaseId) = @_;

    die "You must supply a potential databaseId to nameIsDatabaseId" if !defined($databaseId);

    return exists ($self->{$kIdToStandardName}{$databaseId});

}

############################################################################
sub nameIsAnnotated{
############################################################################
# This method returns a boolean to indicate whether the supplied name has any 
# annotations, either when considered as a databaseId, a standardName, or
# an alias.  If an aspect is also supplied, then it indicates whether that
# name has any annotations in that aspect only.
#
# Usage :
#
# if ($annotationParser->nameIsAnnotated(name => $name)){
#
#       # blah
#
# }
#
# or:
#
# if ($annotationParser->nameIsAnnotated(name   => $name,
#                                        aspect => $aspect)){
#
#       # blah
#
# }

    my ($self, %args) = @_;

    my $name = $args{'name'} || die "You must supply a name to nameIsAnnotated";
    
    my $aspect = $args{'aspect'};

    my $isAnnotated = 0;

    if (!defined ($aspect)){ # if there's no aspect

	$isAnnotated = (exists ($self->{$kNameToIdMap}{$name}) || exists ($self->{$kAmbiguousNames}{$name}));
	
    }else{

	if ($self->nameIsDatabaseId($name) && @{$self->goidsByDatabaseId(databaseId => $name,
									 aspect     => $aspect)}){

	    $isAnnotated = 1;

	}elsif ($self->nameIsStandardName($name) && @{$self->goIdsByStandardName(standardName => $name,
										 aspect       => $aspect)}){

	    $isAnnotated = 1;

	}elsif (!$self->nameIsAmbiguous($name)){

	    my $goidsRef = $self->goIdsByName(name   => $name,
					      aspect => $aspect);

	    if (defined $goidsRef && @{$goidsRef}){

		$isAnnotated = 1;

	    }

	}else { # MUST be an ambiguous name, that's not used as a standard name
	
	    foreach my $databaseId ($self->databaseIdsForAmbiguousName($name)){

		if (@{$self->goidsByDatabaseId(databaseId => $name,
					       aspect     => $aspect)}){

		    $isAnnotated = 1;
		    last; # as soon as we know, we can finish

		}

	    }

	}

    }

    return $isAnnotated;

}

# other public methods

############################################################################
sub databaseName{
############################################################################
# This method returns the name of the annotating authority from the
# file that was supplied to the constructor.
#
# Usage :
#
#    my $databaseName = $annotationParser->databaseName();

    my $self = shift;

    return $self->{$kDatabaseName};

}

############################################################################
sub numAnnotatedGenes{
############################################################################
# This method returns the number of entities in the annotation file
# that have annotations in the supplied aspect.  If no aspect is
# provided, then it will return the number of genes with an annotation
# in at least one aspect of GO.
#
# Usage:
#
#    my $numAnnotatedGenes = $annotationParser->numAnnotatedGenes();
#
#    my $numAnnotatedGenes = $annotationParser->numAnnotatedGenes($aspect);
#

    my ($self, $aspect) = @_;

    if (defined ($aspect)){

	return $self->{$kNumAnnotatedGenes}{$aspect};

    }else{

	return $self->{$kTotalNumAnnotatedGenes};

    }

}

############################################################################
sub allDatabaseIds{
############################################################################
# This public method returns an array of all the database identifiers
#
# Usage:
#
#    my @databaseIds = $annotationParser->allDatabaseIds();
#

    my $self = shift;

    return keys (%{$self->{$kIdToStandardName}});

}

############################################################################
sub allStandardNames{
############################################################################
# This public method returns an array of all standard names.
#
# Usage:
#
#    my @standardNames = $annotationParser->allStandardNames();

    my $self = shift;

    return keys(%{$self->{$kStandardNameToId}});

}

# methods to do with files

############################################################################
sub file{
############################################################################
# This method returns the name of the file that was used to instantiate the
# object.
#
# Usage:
#
#    my $file = $annotationParser->file;

    return $_[0]->{$kFileName};

}

############################################################################
sub serializeToDisk{
############################################################################
# This public method saves the current state of the Annotation Parser
# Object to a file, using the Storable package.  The data are saved in
# network order for portability, just in case.  The name of the object
# file is returned.  By default, the name of the original file will be used
# to make the name of the object file (including the full path from where the
# file came), or the client can instead supply their own filename.
#
# Usage:
#
#    my $fileName = $annotationParser->serializeToDisk;
#    my $fileName = $annotationParser->serializeToDisk(filename => $filename);

    my ($self, %args) = @_;

    my $fileName;

    if (exists ($args{'filename'})){ # they supply their own filename

	$fileName = $args{'filename'};

    }else{ # we build a name from the file used to instantiate ourselves

	$fileName = $self->file;
	
	if ($fileName !~ /\.obj$/){ # if we weren't instantiated from an object
	    
	    $fileName .= ".obj"; # add a .obj suffix to the name
	    
	}

    }

    nstore ($self, $fileName) || die "$PACKAGE could not serialize itself to $fileName : $!";

    return ($fileName);

}

1; # to keep perl happy

############################################################################
#               MORE P O D   D O C U M E N T A T I O N                     #
############################################################################

=pod

=head1 Class Methods

=head2 Usage

This class method simply prints out a usage statement, along with an
error message, if one was passed in.

Usage :

    GO::AnnotationProvider::AnnotationParser->Usage();

=head2 ValidateFile

This class method reads an annotation file, and returns a reference to
an array of errors that are present within the file.  The errors are
simply strings, each beginning with "Line $lineNo : " where $lineNo is
the number of the line in the file where the error was found.

Usage:

    my $errorsRef = GO::AnnotationProvider::AnnotationParser->ValidateFile(annotationFile => $file);

=head1 Constructor

=head2 new

This is the constructor for an AnnotationParser object.

The constructor expects one of two arguments, either a
'annotationFile' argument, or and 'objectFile' argument.  When
instantiated with an annotationFile argument, it expects it to
correspond to an annotation file created by one of the GO consortium
members, according to their file format.  When instantiated with an
objectFile argument, it expects to open a previously created
annotationParser object that has been serialized to disk (see the
serializeToDisk method).

Usage:

    my $annotationParser = GO::AnnotationProvider::AnnotationParser->new(annotationFile => $file);

    my $annotationParser = GO::AnnotationProvider::AnnotationParser->new(objectFile => $file);

=head1 Public instance methods

=head1 Some methods dealing with ambiguous names

Because there are many names by which an annotated entity may be
referred to, that are non-unique, there exist a set of methods for
determining whether a name is ambiguous, and to what database
identifiers such ambiguous names may refer.

=head2 nameIsAmbiguous

This public method returns a boolean to indicate whether a name is
ambiguous, ie whether the name might map to more than one entity (and
therefore more than one databaseId)

Usage:

    if ($annotationParser->nameIsAmbiguous($name)){

        do something useful....or not....

    }

=head2 databaseIdsForAmbiguousName

This public method returns an array of database identifiers for an
ambiguous name.  If the name is not ambiguous, an empty list will be
returned.

Usage:

    my @databaseIds = $annotationParser->databaseIdsForAmbiguousName($name);

=head2 ambiguousNames

This method returns an array of names, which from the annotation file
have been deemed to be ambiguous.

Usage:

    my @ambiguousNames = $annotationParser->ambiguousNames();

=head1 Methods for retrieving GO annotations for entities

=head2 goIdsByDatabaseId

This public method returns a reference to an array of GOIDs that are
associated with the supplied databaseId for a specific aspect.  If no
annotations are associated with that databaseId in that aspect, then a
reference to an empty array will be returned.  If the databaseId is
not recognized, then undef will be returned.

Usage:

    my $goidsRef = $annotationParser->goIdsByDatabaseId(databaseId => $databaseId,
							aspect     => <P|F|C>);

=head2 goIdsByStandardName

This public method returns a reference to an array of GOIDs that are
associated with the supplied standardName for a specific aspect.  If
no annotations are associated with the entity with that standard name
in that aspect, then a a reference to an empty list will be returned.
If the supplied name is not used as a standard name, then undef will
be returned.

Usage:

    my $goidsRef = $annotationParser->goIdsByStandardName(standardName => $databaseId,
							  aspect       => <P|F|C>);


=head2 goIdsByName

This public method returns a reference to an array of GO IDs that are
associated with the supplied name for a specific aspect.  If there are
no GO associations for the entity corresponding to the supplied name
in the provided aspect, then a reference to an empty list will be
returned.  If the supplied name does not correspond to any entity,
then undef will be returned.  Because the name can be any of the
databaseId, the standard name, or any of the aliases, it is possible
that the name might be ambiguous.  Clients of this object should first
test whether the name they are using is ambiguous, using the
nameIsAmbiguous() method, and handle it accordingly.  If an ambiguous
name is supplied, then it will die.

Usage:

    my $goidsRef = $annotationParser->goIdsByName(name   => $name,
						  aspect => <P|F|C>);

=head1 Methods for mapping different types of name to each other

=head2 standardNameByDatabaseId

This method returns the standard name for a database id.

Usage:

    my $standardName = $annotationParser->standardNameByDatabaseId($databaseId);

=head2 databaseIdByStandardName

This method returns the database id for a standard name.

Usage:

    my $databaseId = $annotationParser->databaseIdByStandardName($standardName);


=head2 databaseIdByName

This method returns the database id for any identifier for a gene (eg
by databaseId itself, by standard name, or by alias).  If the used
name is ambiguous, then the program will die.  Thus clients should
call the nameIsAmbiguous() method, prior to using this method.  If the
name does not map to any databaseId, then undef will be returned.

Usage:

    my $databaseId = $annotationParser->databaseIdByName($name);

=head2 standardNameByName

This public method returns the standard name for the the gene
specified by the given name.  Because a name may be ambiguous, the
nameIsAmbiguous() method should be called first.  If an ambiguous name
is supplied, then it will die with an appropriate error message.  If
the name does not map to a standard name, then undef will be returned.

Usage:

    my $standardName = $annotationParser->standardNameByName($name);

=head1 Other methods relating to names
 
=head2 nameIsStandardName

This method returns a boolean to indicate whether the supplied name is
used as a standard name.

Usage :

    if ($annotationParser->nameIsStandardName($name)){

	# do something

    }

=head2 nameIsDatabaseId

This method returns a boolean to indicate whether the supplied name is
used as a database id.

Usage :

    if ($annotationParser->nameIsDatabaseId($name)){

	# do something

    }

=head2 nameIsAnnotated

This method returns a boolean to indicate whether the supplied name has any 
annotations, either when considered as a databaseId, a standardName, or
an alias.  If an aspect is also supplied, then it indicates whether that
name has any annotations in that aspect only.

Usage :

    if ($annotationParser->nameIsAnnotated(name => $name)){

	# blah

    }

or:

    if ($annotationParser->nameIsAnnotated(name   => $name,
					   aspect => $aspect)){

	# blah

    }

=head1 Other public methods

=head2 databaseName

This method returns the name of the annotating authority from the file
that was supplied to the constructor.

Usage :

    my $databaseName = $annotationParser->databaseName();

=head2 numAnnotatedGenes

This method returns the number of entities in the annotation file that
have annotations in the supplied aspect.  If no aspect is provided,
then it will return the number of genes with an annotation in at least
one aspect of GO.

Usage:

    my $numAnnotatedGenes = $annotationParser->numAnnotatedGenes();

    my $numAnnotatedGenes = $annotationParser->numAnnotatedGenes($aspect);

=head2 allDatabaseIds

This public method returns an array of all the database identifiers

Usage:

    my @databaseIds = $annotationParser->allDatabaseIds();

=head2 allStandardNames

This public method returns an array of all standard names.

Usage:

    my @standardNames = $annotationParser->allStandardNames();

=head1 Methods to do with files

=head2 file

This method returns the name of the file that was used to instantiate
the object.

Usage:

    my $file = $annotationParser->file;

=head2 serializeToDisk

This public method saves the current state of the Annotation Parser
Object to a file, using the Storable package.  The data are saved in
network order for portability, just in case.  The name of the object
file is returned.  By default, the name of the original file will be
used to make the name of the object file (including the full path from
where the file came), or the client can instead supply their own
filename.

Usage:

    my $fileName = $annotationParser->serializeToDisk;

    my $fileName = $annotationParser->serializeToDisk(filename => $filename);

=head1 AUTHORS

Elizabeth Boyle, ell@mit.edu

Gavin Sherlock,  sherlock@genome.stanford.edu

=cut
