#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include "TTree.h"
#include "TObjArray.h"
#include "WCPLEEANA/cuts.h"

// B-03: In get_weight, `int wbin;` is declared uninitialized at cuts.h:209.
// In the equal_binning=false path, if no bin edge satisfies the loop
// condition, wbin stays uninitialized and reweight[wbin] reads garbage.
//
// Fix: initialise wbin to a safe value (e.g. 0) and use a `found` flag to
// skip the reweight application when no bin matched.
// Post-fix: when var falls between max bin edge and max_var with no match,
// addtl_weight should remain unchanged (1.0 for "add_weight").

using RwInner = std::tuple<bool, TString, TString, double, double,
                           bool, bool, bool,
                           std::vector<double>, std::vector<double>>;
using RwInfo  = std::tuple<bool, std::vector<RwInner>>;

static RwInfo make_rw_info(bool apply, bool equal_binning,
                            double min_var, double max_var,
                            bool underflow, bool overflow,
                            std::vector<double> reweight,
                            std::vector<double> bins) {
    RwInner entry(apply,
                  "NCPi0",               // cut_str
                  "truth_energyInside",  // var_str
                  min_var, max_var,
                  underflow, overflow,
                  equal_binning,
                  reweight, bins);
    return RwInfo(true, std::vector<RwInner>{entry});
}

TEST_CASE("B-03: no-match bin leaves addtl_weight unchanged") {
    LEEana::EvalInfo   eval   = {};
    LEEana::PFevalInfo pfeval = {};
    LEEana::KineInfo   kine   = {};
    LEEana::TaggerInfo tagger = {};

    // Setup so get_rw_cut_pass("NCPi0") returns true
    eval.truth_isCC        = 0;
    pfeval.truth_NprimPio  = 1;
    pfeval.truth_NCDelta   = 0;

    // var = truth_energyInside = 5.0, which is > min_var=0 and <= max_var=10
    // but no bin in {2.0, 3.0, 4.0} covers 5.0
    eval.truth_energyInside = 5.0f;
    eval.weight_cv     = 1.0f;
    eval.weight_spline = 1.0f;

    auto rw = make_rw_info(
        /*apply*/         true,
        /*equal_binning*/ false,
        /*min_var*/       0.0,
        /*max_var*/       10.0,
        /*underflow*/     false,
        /*overflow*/      false,
        /*reweight*/      {1.5, 2.5, 3.5},
        /*bins*/          {2.0, 3.0, 4.0}   // 5.0 falls in gap after last edge
    );

    // Post-fix: addtl_weight should remain 1.0 when no bin matches
    // (currently UB — stack-garbage index may cause wrong value or crash)
    double result = LEEana::get_weight("add_weight", eval, pfeval, kine, tagger, rw);

    // Expected post-fix: 1.0 (no reweight applied)
    CHECK(result == doctest::Approx(1.0));
}

TEST_CASE("B-03: matching bin applies correct reweight") {
    LEEana::EvalInfo   eval   = {};
    LEEana::PFevalInfo pfeval = {};
    LEEana::KineInfo   kine   = {};
    LEEana::TaggerInfo tagger = {};

    eval.truth_isCC         = 0;
    pfeval.truth_NprimPio   = 1;
    eval.truth_energyInside = 2.5f;  // falls in bin [2.0, 3.0]
    eval.weight_cv          = 1.0f;
    eval.weight_spline      = 1.0f;

    auto rw = make_rw_info(
        true, false, 0.0, 10.0, false, false,
        {1.5, 2.5, 3.5},
        {2.0, 3.0, 4.0}
    );
    // var=2.5 in bin 0 ([2.0,3.0]), reweight[0]=1.5
    double result = LEEana::get_weight("add_weight", eval, pfeval, kine, tagger, rw);
    CHECK(result == doctest::Approx(1.5));
}
