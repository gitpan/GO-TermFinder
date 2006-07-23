use ExtUtils::MakeMaker;
use File::Find;

use strict;
use warnings;

# not clear how to portably identify what c++ compiler there is, so
# we'll see if g++ exists, and if not die, indicating that they need
# to hand edit the $CC variable.

my $CC = 'g++';

my $found = 0; # indicates whether we have found g++

my @directories = split(":", $ENV{"PATH"});

# they may have some garbage in their PATH, so only keep the real
# directories to search

my @realDirectories;

foreach my $directory (@directories){

    if (-e $directory && -d $directory){

	push (@realDirectories, $directory);

    }

}

find(\&wanted, @realDirectories); 

if (!$found){

    # we didn't find $CC

    die "\n\n********************************************************************

*** Error compiling C++ code ***

Not able to locate $CC.

You will have to edit native/Makefile.PL to indicate the name of the
C++ compiler on your system.

********************************************************************
\n\n";

}

# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.

WriteMakefile(
	      'NAME'		=> 'GO::TermFinder::Native',
	      'LIBS'            => ['-lm -lstdc++'],
	      'CC'              => $CC,
	      'LD'              =>'$(CC)',
	      'OBJECT'          => 'Distributions.o Native_wrap.o',
	      );

sub wanted {

    $found = 1 if $_ eq $CC;

}
