#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

TEST_CASE("infrastructure smoke test") {
    CHECK(1 + 1 == 2);
}
