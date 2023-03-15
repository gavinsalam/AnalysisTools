/// tests the usage and I/O of all the main storage types
/// This piece of code also serves as a basic illustration of how to use 
///
/// Usage:
///   ./test-all-types -out outfile -nev 1e4
///
/// Then look at outfile to see what you've got
#include "AnalysisBase.hh"
#include "GSLRandom.hh"

class AnalysisTest : public AnalysisBase {
public:
  double random_number;
  GSLRandom gsl;
  AnalysisTest(CmdLine & cmdline) : AnalysisBase(&cmdline) {}

  void user_startup() override {
    set_default_binning(0.2, 0.8, 0.1);
  }

  /// returns true if the event was successfully generated
  /// (if false, analysis of the event is skipped, but it still
  /// counts towards the total number of events)
  bool generate_event() override {
    random_number = gsl.uniform(0.0, 2.0);
    return true;
  }

  void user_analyse_event() override {
    double evwgt = 1.0;
    xsections["total"] += evwgt;
    averages["RN1"] += evwgt*random_number;

    hists["RN2"].add_entry(random_number, evwgt);
    hists_err["RN3"].add_entry(random_number, evwgt);

    cumul_hists["RN4"].add_entry(random_number, evwgt);
    cumul_hists_err["RN5"].add_entry(random_number, evwgt);

    norm_hists["RN6"].add_entry(random_number, evwgt);
    norm_hists_err["RN7"].add_entry(random_number, evwgt);

    // in the random number bin, get the average of the random number squared
    avg_hists["RN8"].add_entry(random_number, random_number*random_number, evwgt);

    corr_hists["RN9"].add_entry(random_number, random_number*2, random_number*3, evwgt);
  }

};

int main(int argc, char ** argv) {
  CmdLine cmdline(argc, argv);
  AnalysisTest analysis(cmdline);
  analysis.run();
}