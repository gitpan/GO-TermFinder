
use Test::More;

eval "use Test::Pod::Coverage 1.00";

plan skip_all => "This is not an error, Test::Pod::Coverage 1.00 required for testing POD coverage" if $@;

# We exclude the SWIG-generated Native.pm, since it has no POD.

my @modules = Test::Pod::Coverage::all_modules();
my @coverage_modules = ();

for my $module (@modules) {
    push(@coverage_modules, $module) unless $module =~ /.*::Native/;
}

plan tests => scalar(@coverage_modules);

for my $module (@coverage_modules) {
    pod_coverage_ok($module);
}


