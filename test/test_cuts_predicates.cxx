#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include "TTree.h"
#include "TObjArray.h"
#include "WCPLEEANA/cuts.h"

// helpers: build minimal zero-initialized structs and set only the fields each predicate reads

static LEEana::EvalInfo make_eval() {
    LEEana::EvalInfo e = {};
    e.is_match_found_int = false;
    e.is_file_type = false;
    return e;
}

static LEEana::TaggerInfo make_tagger() {
    // Raw-pointer fields stay null; tested predicates only touch scalars.
    LEEana::TaggerInfo t = {};
    return t;
}

// ── is_preselection ──────────────────────────────────────────────────────────

TEST_CASE("is_preselection: passing event") {
    auto e = make_eval();
    e.match_found      = 1;
    e.stm_eventtype    = 1;
    e.stm_clusterlength = 1.0f;
    // stm_lowenergy, stm_LM, stm_TGM, stm_STM, stm_FullDead default to 0
    CHECK(LEEana::is_preselection(e) == true);
}

TEST_CASE("is_preselection: match_found=0 fails") {
    auto e = make_eval();
    e.match_found      = 0;
    e.stm_eventtype    = 1;
    e.stm_clusterlength = 1.0f;
    CHECK(LEEana::is_preselection(e) == false);
}

TEST_CASE("is_preselection: stm_eventtype=0 fails") {
    auto e = make_eval();
    e.match_found      = 1;
    e.stm_eventtype    = 0;
    e.stm_clusterlength = 1.0f;
    CHECK(LEEana::is_preselection(e) == false);
}

TEST_CASE("is_preselection: stm_lowenergy=1 fails") {
    auto e = make_eval();
    e.match_found      = 1;
    e.stm_eventtype    = 1;
    e.stm_lowenergy    = 1;
    e.stm_clusterlength = 1.0f;
    CHECK(LEEana::is_preselection(e) == false);
}

TEST_CASE("is_preselection: clusterlength=0 fails") {
    auto e = make_eval();
    e.match_found      = 1;
    e.stm_eventtype    = 1;
    e.stm_clusterlength = 0.0f;
    CHECK(LEEana::is_preselection(e) == false);
}

TEST_CASE("is_preselection: match_found_asInt path") {
    auto e = make_eval();
    e.is_match_found_int  = true;
    e.match_found_asInt   = 1;
    e.stm_eventtype       = 1;
    e.stm_clusterlength   = 1.0f;
    CHECK(LEEana::is_preselection(e) == true);

    e.match_found_asInt = 0;
    CHECK(LEEana::is_preselection(e) == false);
}

// ── is_FC ────────────────────────────────────────────────────────────────────

TEST_CASE("is_FC: match_isFC=true") {
    auto e = make_eval();
    e.match_isFC = true;
    CHECK(LEEana::is_FC(e) == true);
}

TEST_CASE("is_FC: match_isFC=false") {
    auto e = make_eval();
    e.match_isFC = false;
    CHECK(LEEana::is_FC(e) == false);
}

// ── is_numuCC ────────────────────────────────────────────────────────────────

TEST_CASE("is_numuCC: passing event") {
    auto t = make_tagger();
    t.numu_cc_flag = 0.0f;
    t.numu_score   = 0.91f;
    CHECK(LEEana::is_numuCC(t) == true);
}

TEST_CASE("is_numuCC: flag<0 fails") {
    auto t = make_tagger();
    t.numu_cc_flag = -1.0f;
    t.numu_score   = 0.95f;
    CHECK(LEEana::is_numuCC(t) == false);
}

TEST_CASE("is_numuCC: score<=0.9 fails") {
    auto t = make_tagger();
    t.numu_cc_flag = 1.0f;
    t.numu_score   = 0.9f;  // equal, not strictly >
    CHECK(LEEana::is_numuCC(t) == false);
}

TEST_CASE("is_numuCC: score=0.85 fails") {
    auto t = make_tagger();
    t.numu_cc_flag = 0.0f;
    t.numu_score   = 0.85f;
    CHECK(LEEana::is_numuCC(t) == false);
}

// ── is_nueCC ─────────────────────────────────────────────────────────────────

TEST_CASE("is_nueCC: passing event") {
    auto t = make_tagger();
    t.numu_cc_flag = 0.0f;
    t.nue_score    = 7.01f;
    CHECK(LEEana::is_nueCC(t) == true);
}

TEST_CASE("is_nueCC: flag<0 fails") {
    auto t = make_tagger();
    t.numu_cc_flag = -1.0f;
    t.nue_score    = 8.0f;
    CHECK(LEEana::is_nueCC(t) == false);
}

TEST_CASE("is_nueCC: score=7.0 fails (not strictly >)") {
    auto t = make_tagger();
    t.numu_cc_flag = 0.0f;
    t.nue_score    = 7.0f;
    CHECK(LEEana::is_nueCC(t) == false);
}

TEST_CASE("is_nueCC: score=6.5 fails") {
    auto t = make_tagger();
    t.numu_cc_flag = 0.0f;
    t.nue_score    = 6.5f;
    CHECK(LEEana::is_nueCC(t) == false);
}

// ── is_NC ────────────────────────────────────────────────────────────────────

TEST_CASE("is_NC: cosmict=0 and score<0 passes") {
    auto t = make_tagger();
    t.cosmict_flag = 0.0f;
    t.numu_score   = -0.1f;
    CHECK(LEEana::is_NC(t) == true);
}

TEST_CASE("is_NC: cosmict=1 fails") {
    auto t = make_tagger();
    t.cosmict_flag = 1.0f;
    t.numu_score   = -0.1f;
    CHECK(LEEana::is_NC(t) == false);
}

TEST_CASE("is_NC: score=0 fails (not < 0)") {
    auto t = make_tagger();
    t.cosmict_flag = 0.0f;
    t.numu_score   = 0.0f;
    CHECK(LEEana::is_NC(t) == false);
}

TEST_CASE("is_NC: score>0 fails") {
    auto t = make_tagger();
    t.cosmict_flag = 0.0f;
    t.numu_score   = 0.5f;
    CHECK(LEEana::is_NC(t) == false);
}
