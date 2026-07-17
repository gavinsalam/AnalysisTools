#include "catch.hpp"
#include "SimpleHist2D.hh"
#include "SimpleHist2DWithError.hh"
using namespace Catch::literals;

using namespace std;
  
//TEST_CASE( "SimpleHist2D", "[SimpleHist2D]" ) {
TEMPLATE_TEST_CASE( "SimpleHist2D", "[SimpleHist2D]", SimpleHist2D, SimpleHist2DWithError ) {
  TestType hist(0.0, 10.0, 2.0, 0.0, 10.0, 2.0);

  SECTION("Size constraints") {
    REQUIRE(hist.size() == 25);
    REQUIRE(hist.outflow_size() == 26);
  }  

  SECTION("Filling and querying the histogram") {
    hist.add_entry(5.0, 5.0);
    hist.add_entry(5.0, 5.0, 0.5);
    hist.add_entry(-1.0, 5.0);
    hist.add_entry(5.0, -1.0);

    //                              outflow, total, bin(2,2)
    std::array<double,3> answers = {    2.0,   3.5,     1.5};
    constexpr unsigned i_outflow = 0, i_total = 1, i_bin22 = 2;

    REQUIRE(hist.outflow() == answers[i_outflow]);
    REQUIRE(hist.outflow_bin() == 25);
    REQUIRE(hist.total_weight() == answers[i_total]);
    REQUIRE(hist(2,2) == answers[i_bin22]);
    REQUIRE(hist(2,3) == 0.0_a);

    TestType hist2 = hist*2.0;
    REQUIRE(hist2(2,2) == answers[i_bin22] * 2.0);
    REQUIRE(hist2.outflow() == answers[i_outflow] * 2.0);
    REQUIRE(hist2.total_weight() == answers[i_total] * 2.0);

    TestType hist3 = hist/2.0;
    REQUIRE(hist3(2,2) == answers[i_bin22] / 2.0);
    REQUIRE(hist3.outflow() == answers[i_outflow] / 2.0);
    REQUIRE(hist3.total_weight() == answers[i_total] / 2.0);

    TestType hist4 = hist + hist2;
    REQUIRE(hist4(2,2) == answers[i_bin22] * 3.0);
    REQUIRE(hist4.outflow() == answers[i_outflow] * 3.0);
    REQUIRE(hist4.total_weight() == answers[i_total] * 3.0);
  }
}