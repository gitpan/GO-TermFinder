use Test::More tests => 20;

# this test script simply tests that none of the scripts or libraries
# show any compile time warnings

use strict;
use warnings;
use diagnostics;

$|=1;

# should replace with file find, but this is a very efficient implementation

my @scripts = `find ./examples -type f | xargs grep -l '^#\!/usr/bin/perl'`;

my @libs    = `find ./blib/ -name '*.pm'`;

chomp @scripts;
chomp @libs;

foreach my $file (@scripts, @libs){

    my @lines = `perl -wc $file 2>&1`;

    like($lines[0],  qr/ OK$/, "$file compiles without warnings");

}
