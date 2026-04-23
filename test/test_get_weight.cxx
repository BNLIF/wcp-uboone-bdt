#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include "TTree.h"
#include "TObjArray.h"
#include "WCPLEEANA/cuts.h"

// EvalInfo weight fields are Float_t (32-bit); promote to double to match
// the arithmetic in get_weight.  Derive expected values from floats.
static const float WCV_F  = 1.1f;
static const float WSPL_F = 1.2f;
static const float WLEE_F = 1.3f;
// cv: addtl(double=1.0) * float * float — float args promoted to double each
static const double CV  = 1.0 * (double)WCV_F * (double)WSPL_F;
static const double SPL = (double)WSPL_F;
static const double LEE = (double)WLEE_F;
// lee_spline cases: Float_t * Float_t = float product, then widened to double
static const double LEE_SPL = (double)(WLEE_F * WSPL_F);

static LEEana::EvalInfo make_eval() {
    LEEana::EvalInfo e = {};
    e.weight_cv     = WCV_F;
    e.weight_spline = WSPL_F;
    e.weight_lee    = WLEE_F;
    e.weight_change = 1.0f;
    e.is_match_found_int = false;
    e.is_file_type       = false;
    return e;
}

// Minimal no-reweighting rw_info: first element false → addtl_weight = 1.0
using RwVec = std::vector<std::tuple<bool,TString,TString,double,double,bool,bool,bool,
                                     std::vector<double>,std::vector<double>>>;
using RwInfo = std::tuple<bool, RwVec>;
static RwInfo make_rw_off() {
    return RwInfo(false, RwVec{});
}

static double gw(const char* name) {
    auto e    = make_eval();
    LEEana::PFevalInfo pf = {};
    LEEana::KineInfo   ki = {};
    LEEana::TaggerInfo ta = {};
    auto rw = make_rw_off();
    return LEEana::get_weight(TString(name), e, pf, ki, ta, rw, false);
}

TEST_CASE("get_weight: cv_spline") {
    CHECK(gw("cv_spline") == doctest::Approx(CV).epsilon(1e-9));
}
TEST_CASE("get_weight: cv_spline_cv_spline") {
    CHECK(gw("cv_spline_cv_spline") == doctest::Approx(CV*CV).epsilon(1e-9));
}
TEST_CASE("get_weight: unity") {
    CHECK(gw("unity") == doctest::Approx(1.0).epsilon(1e-9));
}
TEST_CASE("get_weight: unity_unity") {
    CHECK(gw("unity_unity") == doctest::Approx(1.0).epsilon(1e-9));
}
TEST_CASE("get_weight: lee_cv_spline") {
    CHECK(gw("lee_cv_spline") == doctest::Approx(LEE*CV).epsilon(1e-9));
}
TEST_CASE("get_weight: lee_cv_spline_lee_cv_spline") {
    CHECK(gw("lee_cv_spline_lee_cv_spline") == doctest::Approx(LEE*CV * LEE*CV).epsilon(1e-9));
}
TEST_CASE("get_weight: lee_cv_spline_cv_spline") {
    CHECK(gw("lee_cv_spline_cv_spline") == doctest::Approx(LEE * CV*CV).epsilon(1e-9));
}
TEST_CASE("get_weight: cv_spline_lee_cv_spline (alias)") {
    CHECK(gw("cv_spline_lee_cv_spline") == doctest::Approx(LEE * CV*CV).epsilon(1e-9));
}
TEST_CASE("get_weight: spline") {
    CHECK(gw("spline") == doctest::Approx(SPL).epsilon(1e-9));
}
TEST_CASE("get_weight: spline_spline") {
    CHECK(gw("spline_spline") == doctest::Approx(SPL*SPL).epsilon(1e-9));
}
TEST_CASE("get_weight: lee_spline") {
    CHECK(gw("lee_spline") == doctest::Approx(LEE_SPL).epsilon(1e-9));
}
TEST_CASE("get_weight: lee_spline_lee_spline") {
    CHECK(gw("lee_spline_lee_spline") == doctest::Approx(LEE_SPL * LEE_SPL).epsilon(1e-9));
}
TEST_CASE("get_weight: lee_spline_spline") {
    CHECK(gw("lee_spline_spline") == doctest::Approx(LEE * SPL*SPL).epsilon(1e-9));
}
TEST_CASE("get_weight: spline_lee_spline (alias)") {
    CHECK(gw("spline_lee_spline") == doctest::Approx(LEE * SPL*SPL).epsilon(1e-9));
}
TEST_CASE("get_weight: add_weight") {
    CHECK(gw("add_weight") == doctest::Approx(1.0).epsilon(1e-9));
}
TEST_CASE("get_weight: unknown name returns 1") {
    CHECK(gw("nonsense_xyz") == doctest::Approx(1.0).epsilon(1e-9));
}
