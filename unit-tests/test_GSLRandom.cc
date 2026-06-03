#include "GSLRandom.hh"
#include "catch.hpp"

using namespace Catch::literals;
using namespace std;


TEST_CASE("GSLRandom","[GSLRandom]") {
  GSLRandom r;
  REQUIRE(r.half_granularity()/pow(2.0,-33) == 1.0_a);

  GSLRandom rr(r);
  REQUIRE(rr.half_granularity()/pow(2.0,-33) == 1.0_a);

  //----- check GSLRandom::choose(...) -----
  vector<double> prob_choices = {0.5, 1.0, 2.0}; 
  double total = 3.5;
  vector<double> results(3,0.0);
  int nev = 1000000;
  for (int i=0; i<nev; ++i) {
    int bin = r.choose(prob_choices, total);
    results[bin] += 1.0;
  }
  // check within ~ 10sigma (ignoring O(1) normalisation factors in the sigma)
  REQUIRE_THAT(results[0]/nev, Catch::Matchers::WithinAbs(prob_choices[0]/total, 10.0/sqrt(nev)));
  REQUIRE_THAT(results[1]/nev, Catch::Matchers::WithinAbs(prob_choices[1]/total, 10.0/sqrt(nev)));
  REQUIRE_THAT(results[2]/nev, Catch::Matchers::WithinAbs(prob_choices[2]/total, 10.0/sqrt(nev)));
}