#!/usr/bin/env bash
# B-23/L-27: Verify that all Feldman-Cousins method names are correctly spelled.
# Grep for the old typos in header, source, and apps.
REPO="$(cd "$(dirname "$0")/.." && pwd)"
fail=0

for path in "$REPO/inc" "$REPO/src" "$REPO/apps"; do
    count=$(grep -rE 'Fiedman|Fledman' "$path" 2>/dev/null | grep -v '^\s*//' | wc -l)
    if [ "$count" -gt 0 ]; then
        echo "FAIL B-23: misspelled Feldman-Cousins name still present in $path:"
        grep -rnE 'Fiedman|Fledman' "$path" 2>/dev/null | grep -v '^\s*//'
        fail=1
    fi
done

if [ $fail -eq 0 ]; then
    echo "PASS B-23: no Fiedman/Fledman typos in inc/, src/, apps/"
fi

exit $fail
