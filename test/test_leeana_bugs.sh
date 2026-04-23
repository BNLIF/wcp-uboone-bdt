#!/usr/bin/env bash
# Static-analysis tests for L-02, L-03, L-04, L-05.
fail=0

check_absent() {
    local id="$1" file="$2" pattern="$3"
    local count
    count=$(grep -cF "$pattern" "$file" 2>/dev/null; true)
    count=${count:-0}
    if [ "${count:-0}" -gt 0 ] 2>/dev/null; then
        echo "FAIL $id ($file): bad pattern still present: $pattern"
        fail=1
    else
        echo "PASS $id"
    fi
}

check_present() {
    local id="$1" file="$2" pattern="$3"
    local count
    count=$(grep -cF "$pattern" "$file" 2>/dev/null; true)
    count=${count:-0}
    if [ "${count:-0}" -eq 0 ] 2>/dev/null; then
        echo "FAIL $id ($file): required pattern absent: $pattern"
        fail=1
    else
        echo "PASS $id"
    fi
}

# L-02: pot_counting.cxx must not dereference end() when runNo is out of range.
# Old bad pattern: the guard only printed when (it==end && runNo in [4000,50000]),
# so the else-branch ran (dereferencing end()) for any runNo outside that range.
check_absent "L-02-bnb"  "$PWD/apps/pot_counting.cxx" \
    'it == map_bnb_infos.end() && pot.runNo >=4000'
check_absent "L-02-ext"  "$PWD/apps/pot_counting.cxx" \
    'it == map_extbnb_infos.end()  && pot.runNo >=4000'

# L-03: plot_hist.cxx must not use sss(2, sss.Length()-2) (drops last 2 chars).
check_absent "L-03" "$PWD/apps/plot_hist.cxx" 'sss(2, sss.Length()-2)'

# L-04: det_cov_matrix.cxx must use correct spelling "bootstrapping".
check_absent "L-04-typo"   "$PWD/apps/det_cov_matrix.cxx" 'cov_mat_boostrapping_'
check_present "L-04-fixed" "$PWD/apps/det_cov_matrix.cxx" 'cov_mat_bootstrapping_'

# L-05: prune_weightsep24_trees.cxx must guard mcweight->at() with a find() check.
check_present "L-05" "$PWD/apps/prune_weightsep24_trees.cxx" \
    'mcweight->find(knob) == mcweight->end()'

exit $fail
