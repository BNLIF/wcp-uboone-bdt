#!/usr/bin/env bash
# Static-analysis tests for Wave-2 logic/safety fixes.
# Each check uses the same check_absent/check_present pattern as test_leeana_bugs.sh.
REPO="$(cd "$(dirname "$0")/.." && pwd)"
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

# B-06: bayes.cxx destructor f_conv_vec deletion loop must NOT be commented out.
check_absent  "B-06" "$REPO/src/bayes.cxx" \
    '// for (auto it = f_conv_vec.begin(); it!= f_conv_vec.end(); it++){'

# B-08: GPRegressor.cxx must not use raw 'new double[n]' for the fitter parameter array.
check_absent  "B-08" "$REPO/src/GPRegressor.cxx" 'new double[n]'

# B-09: WienerSVD.cxx must explicitly call Decompose() and check its return.
check_present "B-09-decV" "$REPO/src/WienerSVD.cxx" 'decV.Decompose()'
check_present "B-09-udv"  "$REPO/src/WienerSVD.cxx" 'udv.Decompose()'

# B-10/L-11: merge_hist.cxx must close TFile before moving to the next input file.
check_present "B-10" "$REPO/apps/merge_hist.cxx" 'temp_file->Close()'

# B-11: kine.h clear_kine_info must null-check pointers before ->clear().
check_present "B-11" "$REPO/inc/WCPLEEANA/kine.h" \
    'if (tagger_info.kine_energy_particle) tagger_info.kine_energy_particle->clear()'

# B-13: bdt_convert.cxx must not use while(!infile.eof()) anti-pattern.
check_absent  "B-13" "$REPO/apps/bdt_convert.cxx" 'while(!infile.eof())'

# L-09: convert_histo.pl must use string comparison ne, not numeric !=, for "#file" header.
check_present "L-09" "$REPO/LEEana/convert_histo.pl" 'ne "#file"'

# L-14: run_gof.pl must not background stat_pred_cov_matrix before merge_hist needs it.
check_absent  "L-14" "$REPO/LEEana/run_gof.pl" 'stat_pred_cov_matrix -r0 &'

# L-15: pot_counting.cxx must guard argv[2] access with argc > 2.
check_present "L-15" "$REPO/apps/pot_counting.cxx" '(argc > 2)'

# L-18: check_xf.pl must pass $temp[3] (file path column), not $temp[2] (weight index).
check_present "L-18" "$REPO/LEEana/check_xf.pl" 'check_xf_weight_xs $temp[3]'

exit $fail
