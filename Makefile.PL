use ExtUtils::MakeMaker;

# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.

WriteMakefile(
	      'NAME'		=> 'GO-TermFinder',
	      'VERSION_FROM'	=> 'VERSION', # finds $VERSION
	      'PREREQ_PM'	=> { Storable => 1.014 }, # e.g., Module::Name => 1.1
	      ($] >= 5.005 ?    ## Add these new keywords supported since 5.005
	       (ABSTRACT        => 'A collection of modules for working with gene ontologies',
		AUTHOR          => 'Gavin Sherlock <sherlock@genome.stanford.edu>') : ()),
	      );
