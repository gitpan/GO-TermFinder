use ExtUtils::MakeMaker;

# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.

WriteMakefile(
	      'NAME'		=> 'GO-TermFinder',
	      'VERSION_FROM'	=> 'VERSION', # finds $VERSION
	      'PREREQ_PM'	=> { Storable   => '1.0.1.12',
				     GD         => 0, # comment this out if you don't want to use GO::View 
				     GraphViz   => 0, # comment this out if you don't want to use GO::View 
				     CGI        => 0,  # needed for GO::TermFinderReport::Html
				     Test::More => 0
				 },
	      ($] >= 5.005 ?    ## Add these new keywords supported since 5.005
	       (ABSTRACT        => 'A collection of modules for working with gene ontologies',
		AUTHOR          => 'Gavin Sherlock <sherlock@genome.stanford.edu>') : ()),
	      );
