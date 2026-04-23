#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include "TTree.h"
#include "TObjArray.h"
#include "WCPLEEANA/cuts.h"

// is_nueCC: tagger.numu_cc_flag >= 0 && tagger.nue_score > 7.0
// is_FC:    eval.match_isFC
// is_numuCC: tagger.numu_cc_flag >= 0 && tagger.numu_score > 0.9
// flag_truth_inside: vtxX in (-1,254.3], vtxY in (-115,117], vtxZ in (0.6,1036.4]

static LEEana::EvalInfo make_eval_base() {
    LEEana::EvalInfo e = {};
    e.is_match_found_int = false;
    e.is_file_type = false;
    return e;
}

static LEEana::TaggerInfo make_tagger_nueCC() {
    LEEana::TaggerInfo t = {};
    t.numu_cc_flag = 0;
    t.nue_score    = 8.0f;
    return t;
}

static LEEana::TaggerInfo make_tagger_plain() {
    LEEana::TaggerInfo t = {};
    t.numu_cc_flag = 0;
    t.nue_score    = 3.0f;  // not nueCC
    return t;
}

static void set_vtx_inside(LEEana::EvalInfo& e) {
    e.truth_vtxX = 50.0f;
    e.truth_vtxY =  0.0f;
    e.truth_vtxZ = 100.0f;
}

static void set_vtx_outside(LEEana::EvalInfo& e) {
    e.truth_vtxX = -100.0f;  // out of range
}

// KineInfo has raw vector pointers that crash on dereference if null.
// Provide valid empty vectors so is_0p/is_1p/is_0pi loops run safely.
static LEEana::KineInfo make_kine() {
    static std::vector<float> fv;
    static std::vector<int>   iv;
    LEEana::KineInfo ki = {};
    ki.kine_energy_particle = &fv;
    ki.kine_particle_type   = &iv;
    ki.kine_energy_included = &iv;
    ki.kine_energy_info     = &iv;
    return ki;
}

static bool gcp(const char* ch, const char* add_cut, bool flag_data,
                LEEana::EvalInfo& e, LEEana::TaggerInfo& ta) {
    LEEana::PFevalInfo pf = {};
    LEEana::KineInfo   ki = make_kine();
    return LEEana::get_cut_pass(TString(ch), TString(add_cut), flag_data,
                                e, pf, ta, ki);
}

// ── nueCC_FC_bnb : flag_nueCC && flag_FC ─────────────────────────────────────

TEST_CASE("nueCC_FC_bnb: nueCC+FC passes (data)") {
    auto e  = make_eval_base();
    auto ta = make_tagger_nueCC();
    e.match_isFC = true;
    CHECK(gcp("nueCC_FC_bnb", "all", true, e, ta) == true);
}

TEST_CASE("nueCC_FC_bnb: no nueCC fails (data)") {
    auto e  = make_eval_base();
    auto ta = make_tagger_plain();  // nue_score=3 → not nueCC
    e.match_isFC = true;
    CHECK(gcp("nueCC_FC_bnb", "all", true, e, ta) == false);
}

TEST_CASE("nueCC_FC_bnb: no FC fails (data)") {
    auto e  = make_eval_base();
    auto ta = make_tagger_nueCC();
    e.match_isFC = false;
    CHECK(gcp("nueCC_FC_bnb", "all", true, e, ta) == false);
}

// ── LEE_FC_nueoverlay : flag_nueCC && flag_FC && flag_truth_inside ────────────

TEST_CASE("LEE_FC_nueoverlay: nueCC+FC+inside passes") {
    auto e  = make_eval_base();
    auto ta = make_tagger_nueCC();
    e.match_isFC = true;
    set_vtx_inside(e);
    CHECK(gcp("LEE_FC_nueoverlay", "all", false, e, ta) == true);
}

TEST_CASE("LEE_FC_nueoverlay: truth vertex outside fails") {
    auto e  = make_eval_base();
    auto ta = make_tagger_nueCC();
    e.match_isFC = true;
    set_vtx_outside(e);
    CHECK(gcp("LEE_FC_nueoverlay", "all", false, e, ta) == false);
}

// ── nueCC_PC_nueoverlay : flag_nueCC && !flag_FC && flag_truth_inside ─────────

TEST_CASE("nueCC_PC_nueoverlay: nueCC+notFC+inside passes") {
    auto e  = make_eval_base();
    auto ta = make_tagger_nueCC();
    e.match_isFC = false;  // PC = not fully contained
    set_vtx_inside(e);
    CHECK(gcp("nueCC_PC_nueoverlay", "all", false, e, ta) == true);
}

TEST_CASE("nueCC_PC_nueoverlay: FC event fails PC channel") {
    auto e  = make_eval_base();
    auto ta = make_tagger_nueCC();
    e.match_isFC = true;  // FC → !flag_FC is false
    set_vtx_inside(e);
    CHECK(gcp("nueCC_PC_nueoverlay", "all", false, e, ta) == false);
}

// ── BG_nueCC_FC_overlay : nueCC && FC && !(truth nueCC in-FV) ────────────────
// The veto fires when truth_isCC==1 && |nuPdg|==12 && inside → so non-nueCC MC passes

TEST_CASE("BG_nueCC_FC_overlay: non-nueCC-MC passes") {
    auto e  = make_eval_base();
    auto ta = make_tagger_nueCC();
    e.match_isFC  = true;
    e.truth_isCC  = 0;  // NC → veto does NOT fire
    set_vtx_inside(e);
    CHECK(gcp("BG_nueCC_FC_overlay", "all", false, e, ta) == true);
}

TEST_CASE("BG_nueCC_FC_overlay: true nueCC-in-FV vetoed") {
    auto e  = make_eval_base();
    auto ta = make_tagger_nueCC();
    e.match_isFC     = true;
    e.truth_isCC     = 1;
    e.truth_nuPdg    = 12;  // |nuPdg|==12
    set_vtx_inside(e);       // flag_truth_inside=true → veto fires
    CHECK(gcp("BG_nueCC_FC_overlay", "all", false, e, ta) == false);
}

// ── add_cut tokenizer path: "RnumuCCinFV" looked up in map_cuts_flag ─────────

TEST_CASE("BG_nueCC_FC_overlay: add_cut=RnumuCCinFV filters when false") {
    // map_cuts_flag["RnumuCCinFV"] requires nuPdg==14; setting nuPdg=12 → false
    // flag_add becomes false → get_cut_pass returns false before reaching channel code
    auto e  = make_eval_base();
    auto ta = make_tagger_nueCC();
    e.match_isFC  = true;
    e.truth_isCC  = 0;
    e.truth_nuPdg = 12;  // → RnumuCCinFV = false
    set_vtx_inside(e);
    CHECK(gcp("BG_nueCC_FC_overlay", "RnumuCCinFV", false, e, ta) == false);
}
