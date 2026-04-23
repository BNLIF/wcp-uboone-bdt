#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"
#include "WCPLEEANA/GPKernel.h"
#include <cmath>

// B-02: RBFKernel::Mag applies log transforms to p1x/p2x in-place BEFORE
// the zero-length-scale guard, so two identical points with negative
// coordinates produce NaN via log(neg), and NaN != NaN is true in IEEE 754,
// causing the guard to incorrectly return 1e6.
//
// Fix: apply log into a separate array so the guard still sees original coords.

TEST_CASE("B-02: identical points with zero length-scale give Mag=0") {
    // dim 0: length scale = 0, log enabled
    // dims 1-4: length scale = 1, log disabled
    // coeff (par[5]) = 1
    RBFKernel kern(
        {0.0, 1.0, 1.0, 1.0, 1.0, 1.0},
        {true, false, false, false, false}
    );

    // Negative first coordinate triggers log(neg)=NaN pre-fix
    double cx[5] = {-1.0, 1.0, 1.0, 1.0, 1.0};
    GPPoint p1(cx), p2(cx);  // identical

    double mag = kern.Mag(p1, p2);

    // Pre-fix:  log(-1)=NaN, NaN!=NaN → guard fires → returns 1e6  [FAIL]
    // Post-fix: guard uses original -1==-1 → doesn't fire → mag=0  [PASS]
    CHECK(mag == doctest::Approx(0.0));
}

TEST_CASE("B-02: different points with zero length-scale give Mag=1e6") {
    // This case should work correctly both before and after the fix
    // (positive coords, log-monotone, so log(1)!=log(2) iff 1!=2)
    RBFKernel kern(
        {0.0, 1.0, 1.0, 1.0, 1.0, 1.0},
        {true, false, false, false, false}
    );
    double ax[5] = {1.0, 1.0, 1.0, 1.0, 1.0};
    double bx[5] = {2.0, 1.0, 1.0, 1.0, 1.0};
    GPPoint p1(ax), p2(bx);

    double mag = kern.Mag(p1, p2);
    CHECK(mag == doctest::Approx(1e6));
}

TEST_CASE("B-02: zero length-scale, same positive coords give Mag=0") {
    // Should work correctly regardless of fix (positive values, log is monotone)
    RBFKernel kern(
        {0.0, 1.0, 1.0, 1.0, 1.0, 1.0},
        {true, false, false, false, false}
    );
    double cx[5] = {3.0, 2.0, 1.0, 1.0, 1.0};
    GPPoint p1(cx), p2(cx);

    CHECK(kern.Mag(p1, p2) == doctest::Approx(0.0));
}
