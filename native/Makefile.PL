
use ExtUtils::MakeMaker;

# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.

WriteMakefile(
	      'NAME'		=> 'GO::TermFinder::Native',
	      'LIBS'            => ['-lm -lstdc++'],
	      'OBJECT'          => 'Distributions.o Native_wrap.o',
	      );
