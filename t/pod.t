use Test::More;

# This test file will test that all of the pod in any files with a .pm
# or a .pl extension in the distribution have syntactically correct
# pod.

# if the Test::Pod modules are not installed, then the scripts are
# skipped.

eval "use Test::Pod 1.00";

plan skip_all => "This is not an error, Test::Pod 1.00 is required for testing POD" if $@;

all_pod_files_ok();
