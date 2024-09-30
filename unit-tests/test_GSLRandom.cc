#include "GSLRandom.hh"
#include "catch.hpp"

using namespace Catch::literals;
using namespace std;

TEST_CASE("GSLRandom","[GSLRandom]") {
  GSLRandom r;
  REQUIRE(r.half_granularity()/pow(2.0,-33) == 1.0_a);

  GSLRandom rr(r);
  REQUIRE(rr.half_granularity()/pow(2.0,-33) == 1.0_a);
}