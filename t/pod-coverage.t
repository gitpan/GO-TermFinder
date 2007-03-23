#!/usr/bin/perl

use strict;
use warnings;
use diagnostics;

use Test::More;

eval {use Test::Pod::Coverage 1.00};

plan skip_all => "This is not an error, Test::Pod::Coverage 1.00 required for testing POD coverage" if $@;

# We exclude the SWIG-generated Native.pm, since it has no POD.

my @modules = Test::Pod::Coverage::all_modules();

plan tests => scalar(@modules) - 1;

foreach my $module (@modules) {

    next if $module =~ /.*::Native/;

    pod_coverage_ok($module);

}



