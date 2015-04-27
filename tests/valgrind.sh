#!/bin/sh
script_dir=`dirname "$0"`
export GLIBCXX_FORCE_NEW
valgrind \
    --leak-check=full \
    --xml=yes \
    --xml-file=$script_dir/valgrind.log \
    --suppressions=$script_dir/valgrind.supp \
    "$@"

