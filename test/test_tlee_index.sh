#!/bin/bash
# B-04: TLee.cxx Exe_Goodness_of_fit has 7 consecutive "if( index==0 )" blocks.
# Only the last one's settings (axis 0-2600 MeV, 26 bins) ever apply;
# the intended per-channel axis configs for indices 1-6 are dead code.
# Fix: change conditions to if(index==1), if(index==2), ..., if(index==6).
# This test fails (exit 1) while the bug is present, passes after fix.

src="$PWD/src/TLee.cxx"
if [ ! -f "$src" ]; then
    echo "SKIP: $src not found"
    exit 0
fi

count=$(grep -c 'if( index==0 )' "$src")
if [ "$count" -ge 7 ]; then
    echo "FAIL B-04: $src has $count identical 'if( index==0 )' blocks (expected <7 after fix)"
    exit 1
fi
echo "PASS B-04: found $count 'if( index==0 )' blocks"
exit 0
