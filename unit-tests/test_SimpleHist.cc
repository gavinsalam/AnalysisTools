#include "catch.hpp"
#include "SimpleHist.hh"
using namespace Catch::literals;

TEST_CASE( "other", "[other]" ) {
  REQUIRE(2==2);
}

  
TEST_CASE( "SimpleHist", "[SimpleHist]" ) {
  SimpleHist hist(0.0, 10.0, 2.0);

  SECTION("Size constraints") {
    REQUIRE(hist.size() == 5);
    REQUIRE(hist.outflow_size() == 6);
    REQUIRE(hist.binsize() == 2.0_a);
  }

  SECTION("Filling and querying the histogram") {
    hist.add_entry(5.0);
    hist.add_entry(9.0, 0.5);
    hist.add_entry(-1.0);
    hist.add_entry(11.0);
    //
    REQUIRE(hist[2] == 1.0_a);
    REQUIRE(hist[4] == 0.5_a);
    REQUIRE(hist.outflow() == 2.0_a);
    REQUIRE(hist.total_weight() == 3.5_a);
    REQUIRE(hist.n_entries()    == 4.0_a);
    REQUIRE(hist.mean()         == 5.571428571428571_a);
  }

  SimpleHist hist2(hist);
  hist2 += (-1.0)*hist;
  SECTION("checking sum/mult") {
    for (unsigned i = 0; i < hist2.outflow_size(); i++ ){
      DYNAMIC_SECTION("sum/mult Loop round " << i) {
        REQUIRE(hist2[i] == Approx(0.0).margin(1e-10));
      }
    }
  }  
}
