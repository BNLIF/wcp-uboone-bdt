#!/bin/bash
# L-01: Perl pipeline drivers fan-out child processes with trailing '&'
# but omit a final 'wait;' statement, creating a race condition where
# downstream steps start reading files still being written.
# Fix: add 'wait;' after each fan-out block in all affected drivers.
# This test fails (exit 1) while any driver has the bug.

leeana="$PWD/LEEana"
if [ ! -d "$leeana" ]; then
    echo "SKIP: LEEana directory not found"
    exit 0
fi

failed=0
for pl in "$leeana"/*.pl; do
    [ -f "$pl" ] || continue
    # Check if file forks with '&' (line ends in & possibly with spaces/comment)
    if grep -qP '&\s*$' "$pl"; then
        # Require 'wait;' to appear somewhere after the last '&'
        if ! grep -q 'wait;' "$pl"; then
            echo "FAIL L-01: $pl forks with '&' but has no 'wait;'"
            failed=1
        fi
    fi
done

if [ "$failed" -eq 0 ]; then
    echo "PASS L-01: all Perl drivers with '&' fan-out have 'wait;'"
    exit 0
fi
exit 1
