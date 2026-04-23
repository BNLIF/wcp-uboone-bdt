#!/usr/bin/env bash
# B-01: universe-loop seeds must not use j*reweight*77777 (correlated sub-lattice).
# The fix seeds gRandom once before the loop using event-specific data.
src="$PWD/src/master_cov_matrix.cxx"
count=$(grep -cF 'SetSeed(j*reweight*77777)' "$src" 2>/dev/null; true)
count=${count:-0}
if [ "${count:-0}" -gt 0 ] 2>/dev/null; then
    echo "FAIL B-01: $count occurrence(s) of correlated per-universe SetSeed(j*reweight*77777) in $src"
    exit 1
fi
echo "PASS B-01: correlated per-universe seed pattern absent"
