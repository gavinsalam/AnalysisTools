#include "catch.hpp"
#include "SimpleHist.hh"
using namespace Catch::literals;

// a dummy test case.
TEST_CASE( "other", "[other]" ) {
  REQUIRE(2==2);
}

  
TEST_CASE( "SimpleHist", "[SimpleHist]" ) {
  bool verbose = false;
  SimpleHist hist(0.0, 10.0, 2.0);

  SECTION("Size constraints") {
    REQUIRE(hist.size() == 5);
    REQUIRE(hist.outflow_size() == 7);
    REQUIRE(hist.binsize() == 2.0_a);
  }

  SECTION("Filling and querying the histogram") {
    hist.add_entry(5.0);
    hist.add_entry(9.0, 0.5);
    hist.add_entry(-1.0);
    hist.add_entry(-2.0);
    hist.add_entry(11.0);
    // N.N_a should be the same as Approx(N.N)
    REQUIRE(hist[2] == 1.0_a);
    REQUIRE(hist[4] == 0.5_a);
    REQUIRE(hist.underflow() == 2.0_a);
    REQUIRE(hist.overflow() == 1.0_a);
    REQUIRE(hist.outflow() == 3.0_a);
    REQUIRE(hist.total_weight() == 4.5_a);
    REQUIRE(hist.n_entries()    == 5.0_a);
    REQUIRE(hist.mean()         == 3.888888888888888888888888_a);

    if (verbose) {
      hist.output(std::cout) << std::endl; 
      hist.output_diff(std::cout) << std::endl;
      hist.output_cumul(std::cout) << std::endl;
    }
  }

  SimpleHist hist2(hist);
  hist2 += (-1.0)*hist;
  SECTION("checking sum/mult (by requirement of getting zero)") {
    for (unsigned i = 0; i < hist2.outflow_size(); i++ ){
      //DYNAMIC_SECTION("sum/mult Loop round " << i) {
      REQUIRE(hist2[i] == Approx(0.0).margin(1e-10));
        //}
    }
  }

  SECTION("checking histograms with negative range") {
    SimpleHist hist3(10.0, 0.0, 2.0);
    hist3.add_entry(9.0);
    hist3.add_entry(12.0);
    REQUIRE(hist3[0] == 1.0);
    REQUIRE(hist3.underflow() == 1.0);
    REQUIRE(hist3.overflow()  == 0.0);
  }

}
