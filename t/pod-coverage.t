use Test::More;

eval "use Test::Pod::Coverage 1.00";

plan skip_all => "This is not an error, Test::Pod::Coverage 1.00 required for testing POD coverage" if $@;

all_pod_coverage_ok();
