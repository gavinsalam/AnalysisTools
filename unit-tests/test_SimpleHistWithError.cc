#include "catch.hpp"
#include "SimpleHistWithError.hh"
using namespace Catch::literals;

using namespace std;
  
TEST_CASE( "SimpleHistWithError", "[SimpleHistWithError]" ) {
  bool verbose = false;
  SimpleHistWithError hist(0.0, 10.0, 2.0);

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
    REQUIRE(hist.sumsqr(4)      == 0.25_a);
    REQUIRE(hist.sumsqr(2)      == 1.0_a);
    REQUIRE(hist.sumsqr(hist.underflow_bin()) == 2.0_a);
    REQUIRE(hist.sumsqr(hist.overflow_bin())  == 1.0_a);
    REQUIRE(hist.error(2)       == 0.8944271909999159_a);
    REQUIRE(hist.mean()         == 3.888888888888888888888888_a);


    auto hist05 = hist*0.5;
    REQUIRE(hist05.error(2)*2.0 == 0.8944271909999159_a);
    
    if (verbose) {
      for (unsigned i = 0; i < hist.outflow_size(); i++) {
        cout << i << " " << hist[i] << " " << hist.sumsqr(i) << endl;
      }
      hist.output(std::cout) << std::endl; 
      hist.output_diff(std::cout) << std::endl;
      hist.output_cumul(std::cout) << std::endl;
    }
  }


}
