#include "catch.hpp"
#include "AnalysisBase.hh"

using namespace Catch::literals;
using namespace std;

class Analysis : public AnalysisBase {
public:
  Analysis(CmdLine * cmdline_in) : AnalysisBase(cmdline_in) {}

  Binning binning{0.0,2.0,1.0};

  void user_startup() override {
    avg_hists["<y>_v_x"].declare(binning);    

  }

  bool generate_event() override {
    _x = 0.5; 
    _y = 3.0; 
    return true;
  }

  void user_analyse_event() override {
    hists["x"].declare_once(0.0, 2.0, 1.0).add_entry(_x);   
    // this should not change the binning, because it was already declared above 
    avg_hists["<y>_v_x"].declare_once(0.0, 3.0, 1.0).add_entry(_x, _y);    
    avg_hists["<y>_v_x"].declare_once(binning).add_entry(_x, _y);    
  }

private:
  double _x;
  double _y;
};

TEST_CASE( "Analysis", "[Analysis]" ) {
  //bool verbose = false;
  //CmdLine cmdline({"test_analysis", "-nev", "3", "-o", "/tmp/a"});
  CmdLine cmdline({"test_analysis", "-nev", "3", "-o", "/dev/null"});
  Analysis analysis(&cmdline);
  analysis.run();

//  SECTION("Size constraints") {
//    REQUIRE(analysis.hists["x"].size() == 3);
//    REQUIRE(analysis.hists["x"].outflow_size() == 5);
//    REQUIRE(analysis.hists["x"].binsize() == 1.0_a);
//  }
//
//  SECTION("Filling and querying the histogram") {
//    analysis.generate_event();
//    analysis.analyse_event();
//    REQUIRE(analysis.hists["x"][0] == 0.5_a);
//    REQUIRE(analysis.hists["x"].underflow() == 0.0_a);
//    REQUIRE(analysis.hists["x"].overflow() == 0.0_a);
//    REQUIRE(analysis.hists["x"].outflow() == 0.0_a);
//    REQUIRE(analysis.hists["x"].total_weight() == 0.5_a);
//    REQUIRE(analysis.hists["x"].n_entries()    == 1.0_a);
//    REQUIRE(analysis.hists["x"].sumsqr(0)      == 0.25_a);
//    REQUIRE(analysis.hists["x"].sumsqr(1)      == 0.0_a);
//    REQUIRE(analysis.hists["x"].sumsqr(2)      == 0.0_a);
//    REQUIRE(analysis.hists["x"].sumsqr(analysis.hists["x"].underflow_bin()) == 0.0_a);
//    REQUIRE(analysis.hists["x"].sumsqr(analysis.hists["x"].overflow_bin())  == 0.0_a);
//    REQUIRE(analysis.hists["x"].error(0)       == 0.5_a);
//    REQUIRE(analysis.hists["x"].mean()         == 0.5_a);
//  }
}
