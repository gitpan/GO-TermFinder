package GO::View;

#########################################################################
# Module Name  :  View.pm
#
# Date created :  Oct. 2003
# 
# Cared for by Shuai Weng <shuai@genome.stanford.edu>
#
# You may distribute this module under the same terms as perl itself
#########################################################################

# POD documentation - main docs before the code

=pod

=head1 NAME

GO::View


=head1 DESCRIPTION

This perl module generates a graphic that displays the parent and child 
relationships of a selected GO term. It also provides the visualization 
for the GO::TermFinder perl module created by the Stanford Microarray 
Database (SMD). This module is useful when analyzing experimental or 
computational results that produce a set of gene products that may have 
a common function or process.

=head1 SYNOPSIS

    use GO::View;

    my $goView = 

       GO::View->new(-goid               => $goid,
		     -ontologyProvider   => $ontology,
		     -annotationProvider => $annotation,
		     -termFinder         => \@pvalues,
		     -aspect             => 'P',
		     -configFile         => $confFile,
		     -imageDir           => "/tmp",
		     -imageUrlRoot       => "http://www.ABC.com/tmp",
		     -imageName          => "GOview.88.png",
		     -tree               => 'up',
		     -nodeUrl            => $goUrl,
                     -geneUrl            => $geneUrl,
		     -pvalueCutOff       => '0.01',
		     -imageLabel         => "SGD");
				  

    argument              required             expect data and type
    -------------------------------------------------------------------------
    -goid                 No          A gene ontology ID (GOID).
                                      If nothing is passed in, the module 
                                      will use the top goid of each ontology 
                                      branch (i.e, goid for 
				      molecular_function, biological_process,
				      or cellular_component)

    -ontologyProvider	  Yes         An ontology provider instance.

    -annotationProvider   No          An annotation provider instance. It is
                                      required for creating tree for GO Term
                                      Finder result.
    
    -termFinder           No          An array of hash references returned 
                                      from 'findTerms' method of 
                                      GO::TermFinder module. It is required
                                      for creating tree for GO Term Finder 
                                      result. 

    -aspect               No          <P|C|F>. The aspect of the ontology 
                                      provider. It is required for creating 
                                      tree for GO Term Finder result.
    
    -configFile           Yes         The configuration file for setting the
                                      general variables for the graphic 
                                      display. 
				  
    -imageDir             Yes         The directory for storing the newly 
                                      created image file. It must be 
                                      world (nobody) readable and writable
                                      if you want to display the image to 
                                      the web.
 
    -imageUrlRoot         No          The url root for the -imageDir. It is
                                      required if you want to display the
                                      image to the web.

    -imageName            No          The image file name. By default, the 
                                      name will be something like 
                                      'GOview.xxxx.png'. The 'xxxx' will be
                                      the process id.  A suffix for the image (.png
                                      or .gif) should not be provided, as that will
                                      be determined at run time, depending on the
                                      capabilities of the GD library.

    -treeType             No          <up|down>. The tree type. 
                                      
                                      1. up   => display the ancestor tree 
                                                 for the given goid.
                                      2. down => display the descendant tree
                                                 for the given goid.
                                      By default, it will display the 
                                      descendant tree.

    -geneUrl              No          The URL for each Gene to link to.
                                      It needs to have the text <REPLACE_THIS> in 
                                      the url which will be substituted 
                                      by the real goid for a node.

    -nodeUrl              No          The url for each GO node to link to.
                                      It needs to have the text <REPLACE_THIS> in 
                                      the url which will be substituted 
                                      by the real goid for a node.

    -pvalueCutOff         No          The p-value cutoff for displaying
                                      the graphic for GO Term Finder. 
                                      The default is 0.01

    -imageLabel           No          The image label which will appear at
                                      the left bottom corner of the map.
    ------------------------------------------------------------------------

    To display the image on the web:

         $goView->showGraph;
    
    To create and return image file name with full path:
    
         my $imageFile = $goView->createImage;



=head1 FEEDBACK

=head2 Reporting Bugs

Bug reports can be submitted via email 

  shuai@genome.stanford.edu

=head1 AUTHOR

Shuai Weng, shuai@genome.stanford.edu

=head1 COPYRIGHT

Copyright (c) 2003 Stanford University. All Rights Reserved.
This module is free software; you can redistribute it and/or 
modify it under the same terms as Perl itself.


=head1 APPENDIX

The rest of the documentation details each of the public methods. 

=cut

use strict;
use GD;
use GraphViz;

use GO::View::GD;

use vars qw ($PACKAGE $VERSION);

$PACKAGE = 'GO::View';
$VERSION = 0.1;

my $kReplacementText = "<REPLACE_THIS>";

#########################################################################

=head1 METHODS

=cut

#########################################################################
sub new {
#########################################################################

=head2 new

 Title   : new
 Function: Initializes the GO::View object. 
         : Recognized named parameters are -goid, -ontologyProvider,
           -annotationProvider, -termFinder, -aspect, -configFile, 
           -imageDir, -imageUrlRoot, -imageName, -treeType, -nodeUrl, 
           -imageLabel
 Returns : a new object
 Args    : named parameters

=cut

#########################################################################
    my ($class, %args) = @_;

    my $self = bless {}, $class;

    $self->_init(%args);

    return $self;

}

#################################################################
sub graph {
#################################################################
=head2 graph

 Title   : graph
 Usage   : my $graph = $goView->graph;
 Function: Gets the newly created Graphviz instance.   
 Returns : a new Graphviz instance.
 
=cut

################################################################
 
    return $_[0]->{GRAPH};

}

################################################################
sub showGraph {
################################################################

=head2 showGraph

 Title   : showGraph
 Usage   : $goView->showGraph;
 Function: Creates the image and print the map image to a file.  
 Returns : the name of the file to which the image was written
 Throws  : Exception if the imageUrlRoot is not passed to the object.  

=cut
 
#########################################################################
    my ($self) = @_;

    if ($self->graph) {

	$self->_createAndShowImage;
    
    }

    return $self->{IMAGE_FILE};

}

#########################################################################
sub imageFile {
#########################################################################

=head2 imageFile

 Title   : imageFile
 Usage   : my $imageFile = $goView->imageFile;
 Function: Gets the newly created image file name (with full path).  
 Returns : image file name.

=cut
    
#########################################################################

    my ($self) = @_;

    return $self->{IMAGE_FILE};


}

################################################################
sub createImage {
################################################################

=head2 createImage

 Title   : createImage
 Usage   : $goView->createImage; 
 Function: Creates the GO tree image file. Calls it only if you 
           want to create the image file only and do not want to
           display the image.  
 Returns : The newly created image file name with full path.

=cut
    
#########################################################################
    
    my ($self) = @_;

    if ($self->graph) {

	$self->{CREATE_IMAGE_ONLY} = 1;

	return $self->_createAndShowImage;
	
    }

}

################################################################
sub imageMap{
################################################################

=head2 imageMap

 Title    : imageMap
 Usage    : my $map = $goView->imageMap;
 Function : returns the text that constitutes an image map for the
            created image.
 Returns  : a string

=cut

#########################################################################

    return $_[0]->{IMAGE_MAP};

}

################################################################
sub _goid {
################################################################
#
# =head2 _goid 
#
#  Title   : _goid
#  Usage   : my $goid = $self->_goid;
#  Function: Gets the goid that client interface passed in.
#  Returns : GOID
#  Args    : none    
# 
# =cut
# 
#########################################################################

    my ($self) = @_;

    return $self->{GOID};

}

#########################################################################
sub _ontologyProvider {
#########################################################################
#
# =head2 _ontologyProvider 
#
#  Title   : _ontologyProvider
#  Usage   : my $ontology = $self->_ontologyProvider;
#  Function: Gets the ontology provider instance which is passed in by 
#            client interface.
#  Returns : ontology provider instance 
#  Args    : none    
# 
# =cut
# 
#########################################################################

    my ($self) = @_;

    return $self->{ONTOLOGY};

}

#########################################################################
sub _annotationProvider {
#########################################################################
#
# =head2 _annotationProvider 
#
#  Title   : _annotationProvider
#  Usage   : my $annotation = $self->_annotationProvider;
#  Function: Gets the annotation provider instance which is passed in by 
#            client interface.
#  Returns : annotation provider instance 
#  Args    : none    
# 
# =cut
# 
#########################################################################

    my ($self) = @_;

    return $self->{ANNOTATION};

}

#########################################################################
sub _termFinder {
#########################################################################
#
# =head2 _termFinder 
#
#  Title   : _termFinder
#  Usage   : my $termFinder = $self->_termFinderProvider;
#  Function: Gets the term finder result arrayref which is passed in by 
#            client interface.
#  Returns : term finder result arrayref 
#  Args    : none    
# 
# =cut
# 
#########################################################################

    my ($self) = @_;

    return $self->{TERM_FINDER};

}

#########################################################################
sub _init {
#########################################################################
#
# =head2 _init
#
# Title   : _init
# Usage   : n/a; automatically called by new()
# Function: Initializes the variables required for creating the map. 
# Returns : void 
# Args    : n/a
# Throws  : Exception if ontology provider instance or tmp image 
#           directory for storing the image file are not passed 
#           to this object.
#
# =cut
#
#########################################################################

    my ($self, %args) = @_;

    if (!$args{-ontologyProvider} || !$args{-configFile} || 
	!$args{-imageDir}) {

	die "Can't build a $PACKAGE object without the ontologyProvider, configuration file, and tmp image directory name for storing the image file.";
		
    }

    $self->{ONTOLOGY} = $args{-ontologyProvider};

    if ($args{-goid} && $args{-goid} !~ /^[0-9]+$/ && 
	$args{-goid} !~ /^GO:[0-9]+$/) {

	die "The example goid passed to $PACKAGE is GO:0000166.";

    }
	
    my $goid = $args{-goid};

    if ($goid && $goid !~ /^GO:[0-9]+$/) { 

	$goid = $self->_formatGoid($goid); 

    }

    if (!$goid) {  ### create graph for term finder

	# set top goid for the given ontology (molecular_function, 
	# biological_process, or cellular_component).

	$goid = $self->_initGoidFromOntology; 

    }

    $self->{GOID} = $goid;

    # work out the image name and url - note that they will both
    # receive a suffix later on that indicates the type of image that
    # has been output (png or gif)

    $self->{IMAGE_DIR} = $args{-imageDir};

    if ($self->{IMAGE_DIR} !~ /\/$/) {  $self->{IMAGE_DIR} .= "/"; }

    my $suffix;

    if (GD::Image->can('png')){

	$suffix = 'png';

    }else{

	$suffix = 'gif';

    }

    my $imageName;

    if (exists $args{-imageName} && defined $args{-imageName} && $args{-imageName} ne ""){

	$imageName = $args{-imageName};

    }else{

	my $id = $$;

	# now keep incrementing $id, until the name
	# doesn't clash with an existing file

	while (-e $self->{IMAGE_DIR}."GOview.$id.$suffix"){

	    $id++;
     
	}

	$imageName = "GOview.$id";

    }

    $imageName .= ".$suffix";

    $self->{IMAGE_FILE} = $self->{IMAGE_DIR}.$imageName;

    if ($args{-imageUrlRoot}) {
	
	$self->{IMAGE_URL} = $args{-imageUrlRoot};
 
	if ($self->{IMAGE_URL} !~ /\/$/) { $self->{IMAGE_URL} .= "/"; }

	$self->{IMAGE_URL} .= $imageName;

    }else{

	# if we didn't get a root url, we just assume that the image will
	# be in the same directory as the image

	$self->{IMAGE_URL} = "./".$imageName;

    }

    $self->{TREE_TYPE} = $args{-treeType} || 'down';

    my $count;

    if ($args{-annotationProvider} && $args{-termFinder } && 
	$args{-aspect}) {

	$self->{PVALUE_CUTOFF} = $args{-pvalueCutOff} || 0.01;

	$self->{ANNOTATION} = $args{-annotationProvider};

	$self->{TERM_FINDER} = $args{-termFinder};

	$self->{ASPECT} = $args{-aspect};

	$count = $self->_initPvaluesGeneNamesDescendantGoids;

    }
    elsif ($args{-annotationProvider} || $args{-termFinder} || 
	   $args{-aspect}) {

	die "You have to pass annotation provider and term finder instances and GO aspect ([P|F|C]) to $PACKAGE if you want to display graphic for term finder result.";

    }
  
    $self->{IMAGE_LABEL} = $args{-imageLabel};

    $self->{GENE_URL} = $args{-geneUrl};

    $self->{GO_URL} = $args{-nodeUrl};

    $self->_initVariablesFromConfigFile($args{-configFile});

    # only make the graph if we had at least one node passing our p value cutoff

    $self->_createGraph if ($count > 0);

}

################################################################
sub _createGraph {
################################################################
#
# =head2 _createGraph
#
# Title   : _createGraph
# Usage   : $self->_createGraph; 
#           automatically called by _init()
# Function: To create the GraphViz instance and add each descendant/ 
#           ancestor goid into the Graph tree.
# Returns : newly created GraphViz instance.
#
# =cut
#
#########################################################################

    my ($self) = @_;

    ####### If client does not ask for up (ancestor) tree and this is not
    ####### for the special tree paths client asked for the given goid to
    ####### its specified descendants (i.e. for GO Term Finder), then
    ####### we need to determine how many generations of the descendants
    ####### we can display.
    if ($self->{TREE_TYPE} !~ /up/i && 
	!$self->{DESCENDANT_GOID_ARRAY_REF}) {
	
	$self->_setGeneration($self->_goid);

    }

    ###### If this is for the ancestor tree, we need to determine up to 
    ###### which ancestor we want to display the tree paths. We will 
    ###### display tree path up to $self->{TOP_GOID}
    if ($self->{TREE_TYPE} =~ /up/i) {

	$self->_setTopGoid;

	if (!$self->{TOP_GOID}) { return; }
	    
	$self->{NODE_NUM} = 
	    $self->_descendantCount($self->{TOP_GOID},
				    $self->{GENERATION});

    }
    
    my $goid;

    if ($self->{TOP_GOID}) {

	$goid = $self->{TOP_GOID};

    }
    else {

	$goid = $self->_goid;

    }

    $self->_createGraphObject;

    my %foundNode;
    my %foundEdge;
    
    $self->_addNode($goid);
    $foundNode{$goid}++;

    ##### draw go_path for ancestor goid ($self->{GOID}) 
    ##### to each descendant goid in @{$self->{DESCENDANT_GOID_ARRAY_REF}}.
    ##### an example use is for GO Term Finder.
    if ($self->{DESCENDANT_GOID_ARRAY_REF}) {
	
	my $topAncestorNode = 
	    $self->_ontologyProvider->nodeFromId($self->_goid);

	$self->{TERM_HASH_REF_FOR_GOID}{$self->_goid} 
	    = $topAncestorNode->term;

	my $i;

	foreach my $goid (@{$self->{DESCENDANT_GOID_ARRAY_REF}}) {

	    my $childNode = $self->_ontologyProvider->nodeFromId($goid);

	    my @path = $childNode->pathsToAncestor($topAncestorNode);

	    my $found = $self->_addAncestorPathToGraph($childNode, 
						       \@path, 
						       \%foundNode, 
						       \%foundEdge); 
	    
            if (!$found) {
		
		next;

	    }

	    $i++;

	    if ($self->{GENE_NAME_HASH_REF_FOR_GOID} &&
		$self->{GENE_NAME_HASH_REF_FOR_GOID}->{$goid}) {

		my $loci = $self->{GENE_NAME_HASH_REF_FOR_GOID}->{$goid};

		$loci =~ s/:/ /g;

		$loci = "$i:".$loci;

		$self->_addNode($loci);

		$self->_addEdge($goid, $loci);

	    }

	}

	return;

    }
    
    ##### draw part of the tree, and it can go up and down the tree.

    ##### draw up tree and only show ancestors in the paths from the given
    ##### goid to the ancestor goid $self->{TOP_GOID}
    ##### since there are too many nodes...
    if ($self->{TREE_TYPE} =~ /up/i && 
	$self->{NODE_NUM} > $self->{MAX_NODE}) {
	
	my $childNode = $self->_ontologyProvider->nodeFromId($self->_goid);

	my $topAncestorNode = $self->_ontologyProvider->nodeFromId($goid);

	my @path = $childNode->pathsToAncestor($topAncestorNode);

	$self->_addAncestorPathToGraph($childNode, 
				       \@path, 
				       \%foundNode, 
				       \%foundEdge); 

	return;

    }
		
    ##### draw down tree
    my $node = $self->_ontologyProvider->nodeFromId($goid);

    $self->_addChildOfTheNodeToGraph($node, 
				     \%foundNode,
				     \%foundEdge);

    return;
    
}

################################################################
sub _addChildOfTheNodeToGraph {
################################################################
#
# =head2 _addChildOfTheNodeToGraph
#
# Title   : _addChildOfTheNodeToGraph
# Usage   : $self->_addChildOfTheNodeToGraph($node, 
#                                            $foundNodeHashRef,
#                                            $foundEdgeHashRef);
#           automatically called by _createGraph()
# Function: To add each unique descendant of the given node to the 
#           graph tree.
# Returns : void
#
# =cut
#
#########################################################################

    my ($self, $node, $foundNodeHashRef, $foundEdgeHashRef,
	$generation) = @_;

    if (!$generation) { $generation = 1; }
    
    my @childNodes = $node->childNodes;

    foreach my $childNode ($node->childNodes) {

	my $parentGoid = $node->goid;

	my $childGoid = $childNode->goid;

	if (!$$foundNodeHashRef{$parentGoid}) {

	    $self->_addNode($parentGoid);

	    $$foundNodeHashRef{$parentGoid}++;

	}

	if (!$$foundNodeHashRef{$childGoid}) {

	    $self->_addNode($childGoid);

	    $$foundNodeHashRef{$childGoid}++;

	}
	if (!$$foundEdgeHashRef{$parentGoid."::".$childGoid}) {

	    $self->_addEdge($parentGoid, $childGoid);

	    $$foundEdgeHashRef{$parentGoid."::".$childGoid}++;

        }

	if ($generation < $self->{GENERATION}) {

	    $self->_addChildOfTheNodeToGraph($childNode, 
					     $foundNodeHashRef, 
					     $foundEdgeHashRef,
					     $generation++);
	}

    }

}

################################################################
sub _addAncestorPathToGraph {
################################################################
#
# =head2 _addAncestorPathToGraph
#
# Title   : _addAncestorPathToGraph
# Usage   : $self->_addAncestorPathToGraph($node,
#                                          $ancestorPathArrayRef, 
#                                          $foundNodeHashRef,
#                                          $foundEdgeHashRef);
#           automatically called by _createGraph()
# Function: To add each unique ancestor of the given node to the 
#           graph tree.
# Returns : void
#
# =cut
#
#########################################################################

    my ($self, $childNode, $ancestorPathArrayRef, 
	$foundNodeHashRef, $foundEdgeHashRef) = @_; 

    my $found;

    foreach my $ancestorNodeArrayRef (@$ancestorPathArrayRef) {
	
	
	push(@$ancestorNodeArrayRef, $childNode);

	@$ancestorNodeArrayRef = reverse(@$ancestorNodeArrayRef);

	for (my $i = 0; $i < @$ancestorNodeArrayRef; $i++) {

	    my ($goid1, $goid2);

	    if (defined $$ancestorNodeArrayRef[$i]) {

		$goid1 = $$ancestorNodeArrayRef[$i]->goid;

		$self->{TERM_HASH_REF_FOR_GOID}{$goid1} 
		        = $$ancestorNodeArrayRef[$i]->term;

	    }
	    if (defined $$ancestorNodeArrayRef[$i+1]) {

		$goid2 = $$ancestorNodeArrayRef[$i+1]->goid;

		$self->{TERM_HASH_REF_FOR_GOID}{$goid2} 
		        = $$ancestorNodeArrayRef[$i+1]->term;
		
	    }
    
	    if ($goid1 && !$$foundNodeHashRef{$goid1}) {

		$self->_addNode($goid1);

		$$foundNodeHashRef{$goid1}++;

	    }

	    if ($goid1 && $goid2 && 
		!$$foundEdgeHashRef{$goid2."::".$goid1}) {

		$self->_addEdge($goid2, $goid1);

		$$foundEdgeHashRef{$goid2."::".$goid1}++;

	    }

	}

	$found++;

    }

    return $found;
    
}

################################################################
sub _createAndShowImage {
################################################################
#
# =head2 _createAndShowImage
#
# Title   : _createAndShowImage
# Usage   : $self->_createAndShowImage();
#           automatically called by showGraph() and createImage().
# Function: To create the graph tree image. It will print the image to
#           stdout if it is called by showGraph().
# Returns : returns graphText file if the text format is changed. 
#           returns image file name if called by createImage().
#
# =cut
#
#########################################################################

    my ($self) = @_;

    my ($width, $height);

    my $graphText = $self->graph->as_text; 

    if ($graphText =~ /graph \[bb=\"0,0,([0-9]+),([0-9]+)\" *\]\;/) {
       
	$width = $1*$self->{WIDTH_DISPLAY_RATIO};

	$height = $2*$self->{HEIGHT_DISPLAY_RATIO};

    }
    else {

	$self->graph->as_png($self->{IMAGE_DIR}."goPath.$$.png");

        return $self->{IMAGE_DIR}."goPath.$$.png";

    }

    my @graphLine = split(/\n/, $graphText);

    my $border = 25;
    
    my $mapWidth = $width+2*$border;

    my $mapHeight = $height+2*$border;

    my $keyY;

    if ($self->{PVALUE_HASH_REF_FOR_GOID} || 
	!$self->{GENE_NAME_HASH_REF_FOR_GOID}) {

	$keyY = $mapHeight;

	if ($mapWidth < $self->{MIN_MAP_WIDTH}) { 

	    $mapWidth =$self->{MIN_MAP_WIDTH}; 

	}

	if (!$self->{GENE_NAME_HASH_REF_FOR_GOID}) {

	    $mapHeight += int((length($self->{MAP_NOTE})*6/($mapWidth-100))*15) 
		+ 65;

	}
	else {

	    $mapHeight += 50;

	}

    }

    my $gd = GO::View::GD->new(width=>$mapWidth,
				   height=>$mapHeight);

    $gd->im->rectangle(0, 0, $mapWidth-1, $mapHeight-1,
		       $gd->blue);

    $self->_drawFrame($gd, $mapWidth, $mapHeight);

    my (@nodeLine, @edgeLine);

    my $preLine;

    foreach my $line (@graphLine) {

	if ($line =~ /\\$/) { 

	    $line =~ s/\\$//;

	    $line =~ s/^ *//;

	    $preLine .= $line;

	    next;

	}
	elsif ($preLine && $line =~ /\;$/ && 
	       $line !~ / *node[0-9]/) {

	    $line = $preLine.$line;

	    undef $preLine;

	}
	if ($line =~ / *node[0-9]+ *\[(label=.+)\]\;$/i) {

	    push(@nodeLine, $1);

	}
	elsif ($line 
	       =~ / *node[0-9]+ *-> *node[0-9]+ \[pos=\"e,(.+)\"\]\;$/i) {

	    push(@edgeLine, $1);

	}

    }

    if ($self->{PVALUE_HASH_REF_FOR_GOID} && 
	$height > $self->{MIN_MAP_WIDTH_FOR_ONE_LINE_KEY}) {

	### draw keys on the top of the map
	if ($width < $self->{MIN_MAP_WIDTH_FOR_ONE_LINE_KEY}) {

	    $self->{MOVE_Y} = 15;

	}

	$self->_drawKeys($gd, $mapWidth, 5, 'isTop');

    }
    if ($self->{PVALUE_HASH_REF_FOR_GOID} || 
	!$self->{GENE_NAME_HASH_REF_FOR_GOID}) {

	$self->_drawKeys($gd, $mapWidth, $keyY);

    }

    #### draw edges first, then draw nodes
    foreach my $line (@edgeLine) {

	$self->_drawEdge($gd, $height, $border, $line);

    }

    foreach my $line (@nodeLine) {

	$self->_drawNode($gd, $height, $border, $line);

    }

    my $imageFile = $self->{IMAGE_FILE}; 
    my $imageUrl  = $self->{IMAGE_URL};

    open (OUT, ">".$imageFile) || die "Can't create ".$imageFile.":$!";

    binmode OUT;

    if ($gd->im->can('png')) {

	print OUT $gd->im->png;

    }
    else {

	print OUT $gd->im->gif;

    }
    close OUT;

    if ($self->{CREATE_IMAGE_ONLY}) {
	
	return $imageFile;

    }

    if (!$self->{CREATE_IMAGE_ONLY}) {

	my $map = $gd->imageMap;

	if (defined ($map)){

	    $self->{IMAGE_MAP} = 
		
		"<MAP NAME='goPathImage'>\n".
		$gd->imageMap.
		"</MAP>".
		"<center><img src='$imageUrl' usemap='#goPathImage'></center><p>\n";

	}

    }

}

######################################################################
sub _drawNode {
######################################################################
#
# =head2 _drawNode
#
# Title   : _drawNode
# Usage   : $self->_drawNode($gd, $height, $border, $line);
#           automatically called by _createAndShowImage().
# Function: To draw each node.
# Returns : void
#           
# =cut
#
#########################################################################

    my ($self, $gd, $height, $border, $line) = @_;

    if ($line =~ /^label=\"([^\"]*)\", *pos=\"([0-9]+),([0-9]+)\", *width=\"([^\"]+)\", *height=\"([^\"]+)\"/i) {

	my $label = $1;

	my @label = split(/\\n/, $label);

	my $boxW = $4*60;

	my $boxH = $5*60;
	
	$boxH -= 4*(@label-1);

	if ($self->{PVALUE_HASH_REF_FOR_GOID}) {

	    $boxH -= 10;

	}
	
	if (!$self->{MOVE_Y}) { $self->{MOVE_Y} = 0; }

	my $x1 = $2*$self->{WIDTH_DISPLAY_RATIO}-$boxW/2 + $border;

	my $y1 = $height - $3*$self->{HEIGHT_DISPLAY_RATIO} - $boxH/2 + $border + $self->{MOVE_Y};
		
	my $x2 = $x1 + $boxW;

	my $y2 = $y1 + $boxH;

	my $goid;

	if ($label =~ /(GO:[0-9]+)$/) {

	    $goid = $1;

	}

	if (!$goid) {

	    $boxH = 9*(@label) + 4;

	}

	my $geneNum;
	my $totalGeneNum;
	my $barColor;
	my $outline;
	my $linkUrl;

	if ($self->{PVALUE_HASH_REF_FOR_GOID} && $goid) {

	    $barColor = $self->_getBoxColor($gd, $goid);

	    if (!$self->{CREATE_IMAGE_ONLY}) {

		$linkUrl = $self->{GO_URL};

		$linkUrl =~ s/$kReplacementText/$goid/ if $linkUrl;

	    }

	}elsif ($goid) {

	    my $node = $self->_ontologyProvider->nodeFromId($goid);
	    
	    if ($node && $node->childNodes) {
		
		$barColor = $gd->lightBlue;
		
		if (!$self->{CREATE_IMAGE_ONLY}) {

		    $linkUrl = $self->{GO_URL};
		    
		    $linkUrl =~ s/$kReplacementText/$goid/ if $linkUrl;
		    
		}
		
	    }else {
		
		$barColor = $gd->grey;
		
	    }
	    
        }else { 

	    $barColor = $gd->grey;

	}


	my $onInfoText;

	if (($self->{TOP_GOID} && 
	     $goid && $goid eq $self->{TOP_GOID}) ||
	    (!$self->{TOP_GOID} && $goid &&
	     $goid eq $self->_goid)) {

	    $self->_drawUpArrow($gd, $goid, ($x1+$x2)/2-7, 
				($x1+$x2)/2+7, $y1-15, 10, 
				$linkUrl);

	}
	
	#### draw box

	$gd->drawBar(barColor=>$barColor,
		     numX1=>$x1,
		     numX2=>$x2,
		     numY=>$y1,
		     linkUrl=>$linkUrl,
		     barHeight=>$boxH,
		     outline=>$outline,
		     onInfoText=>$onInfoText);
	
	       
	#### draw go_term 

	my $i = 0;

	foreach my $label (@label) {

	    if (!$label || $label =~ /^GO:/i){ next; }

	    if (!$goid) {

		$label =~ s/[0-9]+://i;

	    }

	    my $nameColor = $gd->black;

	    if (!$goid) {

		$nameColor = $gd->blue;

	    }
	    elsif ($goid eq $self->_goid) {

		$nameColor = $gd->red;

	    }

	    my $startPixel = int(($boxW - length($label)*6)/2);
	    my $numX1 = $x1 + $startPixel;
	    my $numY1 = $y1 + $i*10;

	    if ($goid) {

		$gd->drawName(name=>$label,
			      nameColor=>$nameColor,  
			      numX1=>$numX1,
			      numY=>$numY1);

	    }
	    else {

		# $numX1 -= 10;

		my @gene = split(' ', $label);

		foreach my $gene(@gene) {
		    
		    my $linkUrl;

		    if (!$self->{CREATE_IMAGE_ONLY} && $self->{GENE_URL}) {

			$linkUrl = $self->{GENE_URL};
			
			$linkUrl =~ s/$kReplacementText/$gene/;

		    }
		    $gd->drawName(name=>$gene,
				  nameColor=>$nameColor,
				  linkUrl=>$linkUrl,
				  numX1=>$numX1,
				  numY=>$numY1);

		    $numX1 += (length($gene)+1)*6;

		}

	    }

	    $i++;

	}

	if ($geneNum) {

	    my $label = $geneNum." gene";

	    if ($totalGeneNum != 1) {

		$label .= "s";

	    }

	    my $startPixel = int(($boxW - length($label)*6)/2);

	    my $numX1 = $x1 + $startPixel+2;

	    my $numY1 = $y1 + $i*10+2;

	    $gd->drawName(name=>$label,
			  nameColor=>$gd->maroon,  
			  numX1=>$numX1,
			  numY=>$numY1);

	}

    }

}

######################################################################
sub _drawEdge {
######################################################################
#
# =head2 _drawEdge
#
# Title   : _drawEdge
# Usage   : $self->_drawEdge($gd, $height, $border, $line);
#           automatically called by _createAndShowImage().
# Function: To draw each edge.
# Returns : void
#           
# =cut
#
#########################################################################

    my ($self, $gd, $height, $border, $line) = @_;

    my @point = split(/ /, $line);
    
    shift @point;

    my ($preX, $preY);

    foreach my $point (@point) {

	my ($x, $y) = split(/\,/, $point);

	$x *= $self->{WIDTH_DISPLAY_RATIO};

	$x += $border;

	if (!defined $self->{MOVE_Y}) { $self->{MOVE_Y} = 0; }

	$y = $height - $y*$self->{HEIGHT_DISPLAY_RATIO} + $border +
	     $self->{MOVE_Y} + 5;

	if ($preX && $preY) {

	    $gd->im->line($x, $y, $preX, $preY, $gd->black);

	}

	$preX = $x;

	$preY = $y;

    }

}

#################################################################
sub _drawUpArrow {
#################################################################
#
# =head2 _drawUpArrow
#
# Title   : _drawUpArrow
# Usage   : $self->_drawUpArrow($gd, $goid, $X1, $X2, $Y, $barHeight,
#                               $linkUrl);
#           automatically called by _drawNode().
# Function: To draw an up arrow on the tree map. 
# Returns : void
#           
# =cut
#
#########################################################################

    my ($self, $gd, $goid, $X1, $X2, $Y, $barHeight, 
	$linkUrl) = @_;

    my $node = $self->_ontologyProvider->nodeFromId($goid);

    my $maxGenerationUp = $node->lengthOfLongestPathToRoot;

    if ($maxGenerationUp <= 1) { return; } 

    if (!$self->{CREATE_IMAGE_ONLY} && $linkUrl) {
    
	$linkUrl .= "&tree=up";
    
    }

    $gd->drawBar(barColor=>$gd->blue,
		 numX1=>$X1,
		 numX2=>$X2,
		 numY=>$Y,
		 linkUrl=>$linkUrl,
		 barHeight=>$barHeight,
		 outline=>1,
		 arrow=>'up');

}

######################################################################
sub _drawKeys {
######################################################################
#
# =head2 _drawKeys
#
# Title   : _drawKeys
# Usage   : $self->_drawKeys($gd, $mapWidth, $keyY, $isTop);
#           automatically called by _createAndShowImage().
# Function: To draw the display keys on the top or bottom of the tree map. 
# Returns : void
#           
# =cut
#
#########################################################################

    my ($self, $gd, $mapWidth, $keyY, $isTop) = @_;

    if (!$self->{GENE_NAME_HASH_REF_FOR_GOID}) {

	my $y = $keyY;

	my $boxH = 10;

	my $startX = 50;
	
	$gd->drawBar(barColor=>$gd->lightBlue,
		     numX1=>$startX,
		     numX2=>$startX+20,
		     numY=>$y,
		     barHeight=>$boxH);

	my $numX1 = $startX + 20;

	$gd->drawName(name=>" = GO term with child(ren)",
		      nameColor=>$gd->black,  
		      numX1=>$numX1,
		      numY=>$y-2);

	$y += 15;

	$gd->drawBar(barColor=>$gd->grey,
		     numX1=>$startX,
		     numX2=>$startX+20,
		     numY=>$y,
		     barHeight=>$boxH);

	$gd->drawName(name=>" = GO term with no child(ren)",
		      nameColor=>$gd->black,  
		      numX1=>$numX1,
		      numY=>$y-2);

	my $maxTextLen = int(($mapWidth-2*$startX)/6);

	my $text = $self->_processLabel($self->{MAP_NOTE}, $maxTextLen);

	my @geneNumExample = split(/\12/, $text);
	
	$y += 15;

	foreach my $text (@geneNumExample) {

	    $text =~ s/^ *//;

	    $gd->drawName(name=>$text,
			  nameColor=>$gd->black,  
			  numX1=>$startX,
			  numY=>$y-2);

	    $y += 15;

	}
 
	return;

    }

    my $y1 = $keyY;
    my $boxH = 15;
    my $boxW = 88;
    my $startX = 48 + ($mapWidth - $boxW*6 - 35 - 48)/2;
    
    my $twoLine;
    if ($startX < 48) {
	$startX = 48 + ($mapWidth - $boxW*3 - 15 - 48)/2;
	$twoLine = 1;
	if (!$isTop) { $y1 -= 10; }
    }

    ####### new code

    if ($isTop) {

	#$startX = 48;

    }

    $gd->drawName(name=>'pvalue:',
		  nameColor=>$gd->black,  
		  numX1=>10,
		  numY=>$y1+1);
    
    my $i;
    my $preX2 = $startX;

    foreach my $name ('<=1e-10', '1e-10 to 1e-8', '1e-8 to 1e-6', 
			'1e-6 to 1e-4', '1e-4 to 1e-2', '>0.01') {
	
	$i++;

	my $pvalue = $name;

	$pvalue =~ s/^<=//;

	$pvalue =~ s/^>0.01/1/;

	$pvalue =~ s/^.+ to (.+)$/$1/;

	my $barColor = $self->_color4pvalue($gd, $pvalue);

	if ($i == 4 && $twoLine) {
	    $preX2 = $startX;
	    $y1 += 20;
	}

	my $x1 = $preX2 + 5; 

	my $x2 = $x1 + $boxW;

	$gd->drawBar(barColor=>$barColor,
		     numX1=>$x1,
		     numX2=>$x2,
		     numY=>$y1,
		     barHeight=>$boxH);

	my $numX1 = $x1 + ($boxW-length($name)*6)/2;

        my $numY1 = $y1 + 2;

	$gd->drawName(name=>$name,
		      nameColor=>$gd->black,  
		      numX1=>$numX1,
		      numY=>$numY1);

	$preX2 = $x2;

    }

}

######################################################################
sub _drawFrame {
######################################################################
#
# =head2 _drawFrame
#
# Title   : _drawFrame
# Usage   : $self->_drawFrame($gd, $width, $height);
#           automatically called by _createAndShowImage().
# Function: To draw a frame around the image map with date and label 
#           on the bottom corner. 
# Returns : void
#           
# =cut
#
#########################################################################

    my ($self, $gd, $width, $height) = @_;

    $gd->drawFrameWithLabelAndDate(width=>$width,
				   height=>$height,
				   text=>$self->{IMAGE_LABEL});

}

################################################################
sub _createGraphObject {
################################################################
#
# =head2 _createGraphObject
#
# Title   : _createGraphObject
# Usage   : my $self->_createGraphObject();
#           automatically called by _createGraph().
# Function: Gets newly created empty GraphViz instance. 
# Returns : newly created empty GraphViz instance.
#           
# =cut
#
#########################################################################

    my ($self) = @_;

    $self->{GRAPH} = GraphViz->new(node=>{shape=>'box'});

}

################################################################
sub _addNode {
################################################################
#
# =head2 _addNode 
#
# Title   : _addNode 
# Usage   : $self->_addNode($goid);
#           automatically called by _createGraph().
# Function: Adds node to the GraphViz instance. 
# Returns : void
#           
# =cut
#
#########################################################################

    my ($self, $goid) = @_;

    if ($goid !~ /^GO:[0-9]+$/) {

	my $label = $self->_processLabel($goid, 30);

	$self->graph->add_node($goid,
			 label=>$label);

	return;

    }

    my $label = $self->{TERM_HASH_REF_FOR_GOID}->{$goid};

    if (!$label) {

	my $node = $self->_ontologyProvider->nodeFromId($goid);

	$label = $node->term;

    }

    my $stdGoid;

    if (!$self->{PVALUE_HASH_REF_FOR_GOID}) {

	$stdGoid = $self->_formatGoid($goid);

    }
    else {

	$stdGoid = "GO:".$goid;

    }

    if ($label) {

	$label = $self->_processLabel($label)."\n".$stdGoid;

    }
    else { $label = $stdGoid; }
   
    $self->graph->add_node($goid,
			   label=>$label);
  
    return;

}

################################################################
sub _addEdge {
################################################################
#
# =head2 _addEdge 
#
# Title   : _addEdge 
# Usage   : $self->_addEdge($parentGoid, $childGoid);
#           automatically called by _createGraph().
# Function: Adds edge to the GraphViz instance. 
# Returns : void
#           
# =cut
#
#########################################################################

    my ($self, $parentGoid, $childGoid) = @_;

    $self->graph->add_edge($parentGoid=>$childGoid);

}

########################################################################
sub _descendantCount {
########################################################################
#
# =head2 _descendantCount
#
# Title   : _descendantCount
# Usage   : my $nodeCount = 
#                $self->_descendantCount($goid, $generationDown);
#           automatically called by _createGraph().
# Function: Gets total descendant node number down to a given generation. 
# Returns : The total descendant node number down to a given generation.
#           
# =cut
#
#########################################################################

    my ($self, $goid, $generationDown) = @_;

    my $node = $self->_ontologyProvider->nodeFromId($goid);

    my %descendantCount4generation;

    $self->_descendantCount4generation($node, \%descendantCount4generation); 

    my $nodeNum = 0;
    
    foreach my $generation (1..$generationDown) {

	$nodeNum += $descendantCount4generation{$generation};
	
    }

    return $nodeNum;
    
}

########################################################################
sub _setGeneration {
#######################################################################
#
# =head2 _setGeneration
#
# Title   : _setGeneration
# Usage   : $self->_setGeneration($goid);
#           automatically called by _createGraph().
# Function: Sets the maximum generation number it will display.
# Returns : void
#           
# =cut
#
#########################################################################

    my ($self, $goid) = @_;
 
    my $node = $self->_ontologyProvider->nodeFromId($goid);

    my %descendantCount4generation;

    $self->_descendantCount4generation($node, \%descendantCount4generation); 
    
    my $nodeNum = 0; 
    my $preNodeNum = 0;
    
    foreach my $generation (sort {$a<=>$b} (keys %descendantCount4generation)) {

	$nodeNum += $descendantCount4generation{$generation};

	if ($nodeNum == $self->{MAX_NODE}) { 

	    $self->{GENERATION} = $generation;

	    $self->{NODE_NUM} = $nodeNum;

	    last;

	}

	if ($nodeNum > $self->{MAX_NODE}) {

	    $self->{GENERATION} = $generation-1;

	    $self->{NODE_NUM} = $preNodeNum;

	    last;

	}
	
        $preNodeNum = $nodeNum;

    }

    if (!$self->{GENERATION} || $self->{GENERATION} < 1) { 

	$self->{GENERATION} = 1;

	if (!$node->childNodes) {

	    $self->{TREE_TYPE} = 'up';
 
	}

    }

}

########################################################################
sub _descendantCount4generation {
########################################################################
#
# =head2 _descendantCount4generation
#
# Title   : _descendantCount4generation
# Usage   : $self->_descendantCount4generation($node, 
#                                              $nodeCountHashRef,
#                                              $generation);
#           automatically called by _descendantCount(),
#                                   _setGeneration(), and itself.
# Function: Gets the descebdant count for each generation.
# Returns : void
#           
# =cut
#
#########################################################################

    my ($self, $node, $nodeCountHashRef, $generation) = @_;
    
    if (!$generation) { $generation = 1; }
    
    my @childNodes = $node->childNodes;

    $$nodeCountHashRef{$generation} += scalar(@childNodes);

    foreach my $childNode ($node->childNodes) {

	$self->_descendantCount4generation($childNode, 
					   $nodeCountHashRef, 
					   $generation++);

    }

}

################################################################
sub _setTopGoid {
################################################################
#
# =head2 _setTopGoid
#
# Title   : _setTopGoid
# Usage   : $self->_setTopGoid();
#           automatically called by _createGraph().
# Function: Sets the top ancestor goid. We want to display the 
#           tree up to this goid.  
# Returns : void
#           
# =cut
#
#########################################################################

    my ($self) = @_;

    my $node = $self->_ontologyProvider->nodeFromId($self->_goid);

    my $maxGenerationUp = $node->lengthOfLongestPathToRoot - 1;

    my @pathsToRoot = $node->pathsToRoot;

    my %count4goid;
    my %generation4goid;

    my $pathNum = scalar(@pathsToRoot);

    foreach my $path (@pathsToRoot) {

	my @nodeInPath = reverse(@$path);

	my $generation = 0;

	foreach my $node (@nodeInPath) {

	    $count4goid{$node->goid}++;

	    $generation++;

	    if (!$generation4goid{$node->goid} || 
		$generation4goid{$node->goid} < $generation) {

		$generation4goid{$node->goid} = $generation;

	    }

	    if ( $pathNum == $count4goid{$node->goid} ) {
		### the same goid appears on all the paths.
		
		$self->{TOP_GOID} = $node->goid;
		
		$self->{GENERATION} = $generation4goid{$node->goid};

		last;

	    }
 
	}

    }

}

################################################################
sub _initGoidFromOntology {
################################################################
#
# =head2 _initGoidFromOntology
#
# Title   : _initGoidFromOntology
# Usage   : my $goid = $self->_initGoidFromOntology;
#           automatically called by _init().
# Function: Gets the top goid for the given ontology  
#           (biological_process, molecular_function, or 
#           cellular_component)
# Returns : goid
#           
# =cut
#
#########################################################################

    my ($self) = @_;

    # gene_ontology super node
    my $rootNode = $self->_ontologyProvider->rootNode;

    # node for molecular_function, biological_process, or 
    # cellular_component 
    my ($topNode) = $rootNode->childNodes;

    return $topNode->goid;

}

#######################################################################
sub _initPvaluesGeneNamesDescendantGoids {
#######################################################################
#
# =head2 _initPvaluesGeneNamesDescendantGoids
#
# Title   : _initPvaluesGeneNamesDescendantGoids
# Usage   : $self->_initPvaluesGeneNamesDescendantGoids;
#           automatically called by _init().
# Function: Sets $self->{PVALUE_HASH_REF_FOR_GOID}, 
#           $self->{GENE_NAME_HASH_REF_FOR_GOID},
#           and $self->{DESCENDANT_GOID_ARRAY_REF} 
# Returns : void
#           
# =cut
#
#########################################################################

    my ($self) = @_;
 
    my %foundGoid;
    my %foundGoidGene;

    my @directAnnotatedGoid;
    my %loci4goid;
    my %pvalue4goid;

    my $maxTopNodeToShow = 6;

    my $count = 0;

    foreach my $pvalue (@{$self->_termFinder}){
    
	next if ($pvalue->{CORRECTED_PVALUE} > $self->{PVALUE_CUTOFF});
        
	last if ($count >= $maxTopNodeToShow);
 
	$pvalue4goid{$pvalue->{NODE}->goid} = $pvalue->{CORRECTED_PVALUE};

	my $ancestorNode = 
	    $self->_ontologyProvider->nodeFromId($pvalue->{NODE}->goid);

	foreach my $databaseId (keys %{$pvalue->{ANNOTATED_GENES}}) {

	    my $gene = $pvalue->{ANNOTATED_GENES}->{$databaseId};

	    my $goidArrayRef = 
		$self->_annotationProvider->goIdsByName(name=>$gene,
						        aspect=>$self->{ASPECT});
	
	    foreach my $goid (@$goidArrayRef) {

		my $node = $self->_ontologyProvider->nodeFromId($goid);
	
		if (!$node) { next; }

		if ($ancestorNode->goid ne $goid && 
		     !$ancestorNode->isAnAncestorOf($node)) {

		    next;

	        }
		if (!$foundGoidGene{$goid."::".$gene}) {  

		    $loci4goid{$goid} .= ":".$gene;

		    $foundGoidGene{$goid."::".$gene}++;

		}

		if ($foundGoid{$goid}) { next; }

		push(@directAnnotatedGoid, $goid);
	    
		$foundGoid{$goid}++;
	    }
	}

	$count++;
   
    }

    $self->{DESCENDANT_GOID_ARRAY_REF} = \@directAnnotatedGoid;
    $self->{PVALUE_HASH_REF_FOR_GOID} = \%pvalue4goid;
    $self->{GENE_NAME_HASH_REF_FOR_GOID} = \%loci4goid;

    return $count;

}

##########################################################################
sub _initVariablesFromConfigFile {
##########################################################################
    my ($self, $configFile) = @_;

    open(CONF, "$configFile") || 
	die "Can't open '$configFile' for reading:$!";

    while(<CONF>) {

	chomp;

	# skip comments, blank, and whitespace only lines

	if (/^\#/ || /^\s*$/) { next; }
	
	my ($name, $value) = split(/=/);

	$value =~ s/^ *(.+) *$/$1/;
	
	if ($name =~ /^maxNode/i) {

	    $self->{MAX_NODE} = $value;

	}
	if ($name =~ /^maxNodeNameWidth/i) {

	    $self->{MAX_NODE_NAME_WIDTH} = $value;

	}
	if ($name =~ /^widthDisplayRatio/i) {

	    $self->{WIDTH_DISPLAY_RATIO} = $value;

	}
	if ($name =~ /^heightDisplayRatio/i) {

	    $self->{HEIGHT_DISPLAY_RATIO} = $value;

        }
	if ($name =~ /^minMapWidth/i) {

	    $self->{MIN_MAP_WIDTH} = $value;

	}
	if ($name =~ /^minMapHeight4TopKey/i) {

	    $self->{MIN_MAP_HEIGHT_FOR_TOP_KEY} = $value;

	}
	if ($name =~ /^minMapWidth4OneLineKey/i) {

	    $self->{MIN_MAP_WIDTH_FOR_ONE_LINE_KEY} = $value;

	}
	if ($name =~ /^mapNote/i) {

	    $self->{MAP_NOTE} = $value;
    
	}
	if ($name =~ /^binDir/i) {

	    $ENV{PATH} .= ":".$value;

	}
	if ($name =~ /^libDir/i) {

	    $ENV{LD_LIBRARY_PATH} .= ":".$value;

	}
    
    }
    close(CONF);

}

################################################################
sub _getBoxColor {
################################################################
#
# =head2 _getBoxColor
#
# Title   : _getBoxColor
# Usage   : my $boxColor = $self->_getBoxColor($gd, $goid);
#           automatically called by _drawNode().
# Function: Gets the color for the node box in the display.
# Returns : gd color
#           
# =cut
#
#########################################################################

    my ($self, $gd, $goid) = @_;
    
    if ($self->{PVALUE_HASH_REF_FOR_GOID} && 
	$self->{PVALUE_HASH_REF_FOR_GOID}->{$goid}) {

	return $self->_color4pvalue($gd, 
				    $self->{PVALUE_HASH_REF_FOR_GOID}->{$goid});
     
    }

    return $gd->tan;

}

################################################################
sub _color4pvalue {
################################################################
#
# =head2 _color4pvalue
#
# Title   : _color4pvalue
# Usage   : my $boxColor = $self->_color4pvalue($gd, $pvalue);
#           automatically called by _drawKeys() and _getBoxColor().
# Function: Gets the color for the node box in the display.
# Returns : gd color
#           
# =cut
#
#########################################################################

    my ($self, $gd, $pvalue) = @_;

    if ($pvalue <= 1e-10) {

	return $gd->orange; 

    }
    elsif ($pvalue <= 1e-8) {

	return $gd->yellow; 

    }
    elsif ($pvalue <= 1e-6) {

	return $gd->green4;

    }
    elsif ($pvalue <= 1e-4) {

	return $gd->lightBlue;

    }
    elsif ($pvalue <= 1e-2) {

	return $gd->blue4;

    }
    else {

	return $gd->tan;

    }

}

################################################################
sub _processLabel {
################################################################
#
# =head2 _processLabel
#
# Title   : _processLabel
# Usage   : my $newLabel = $self->_processLabel($label,
#                                               $maxLabelLen);
#           automatically called by _drawKeys() and _addNode().
# Function: Splits the label into multiple lines if the label is too 
#           long. 
# Returns : new label string
#           
# =cut
#
#########################################################################

    my ($self, $label, $maxLabelLen) = @_;
    
    if (!$maxLabelLen) { 

	$maxLabelLen = $self->{MAX_NODE_NAME_WIDTH} || 15;

    }

    my @word = split(/ /, $label);

    undef $label;

    my $line;

    foreach my $word (@word) {

	if ( ($line && length($line) >= $maxLabelLen) ||
	     ($line && 
	      (length($line)+length($word) > $maxLabelLen)) ) {

	    $line =~ s/^ +//;

	    $label .= "\n".$line;

	    undef $line;

	}

	$line .= " ".$word;    

    }
    if ($line) {

	$label .= "\n".$line;

    }

    $label =~ s/^\n//;

    return $label;

}

################################################################
sub _formatGoid {
################################################################
#
# =head2 _formatGoid
#
# Title   : _formatGoid
# Usage   : my $goid = $self->_formatGoid($goid); 
#           automatically called by _init() and _addNode().
# Function: Reformats the goid (plain number) to STD GOID format 
#           (GO:0000388)
# Returns : std GOID
#           
# =cut
#
#########################################################################

    my ($self, $goid) = @_;

    my $len = length($goid);

    for (my $i = 1; $i <= 7 - $len; $i++) {

	$goid = "0".$goid;

    }

    $goid = "GO:".$goid;

    return $goid;

}

#######################################################################
sub DESTROY {
#######################################################################

    # nothing needs to be done

}

#######################################################################
1;
#######################################################################

