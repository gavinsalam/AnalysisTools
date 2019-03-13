#include "AnalysisBase.hh"
#include <fstream>

using namespace std;

template<> double DefaultHist::_lo =   0.0;
template<> double DefaultHist::_hi = 100.0;
template<> double DefaultHist::_bin_size = 1.0;

template<> double DefaultAveragingHist::_lo =   0.0;
template<> double DefaultAveragingHist::_hi = 100.0;
template<> double DefaultAveragingHist::_bin_size = 1.0;

template<> double DefaultCorrelationHist::_lo =   0.0;
template<> double DefaultCorrelationHist::_hi = 100.0;
template<> double DefaultCorrelationHist::_bin_size = 1.0;


//======================================================================
AnalysisBase::AnalysisBase(CmdLine * cmdline_in) : cmdline(cmdline_in) {

  header << cmdline->header() << endl;
  nev = cmdline->value<double>("-nev", 1e2);
  output_filename = cmdline->value<string>("-out");
  
}

//======================================================================
void AnalysisBase::run() {
  pre_startup();
  user_startup();
  post_startup();
  user_post_startup();

  cmdline->assert_all_options_used();

  /// should this go into a separate event_loop routine?
  for (iev = 0; iev < nev; ) {

    bool success = generate_event();
    // what do we do on failure: bail out, or continue round the loop?
    if (! success) break;
    
    total_weight += event_weight();

    user_analyse_event();

    // certain histogram types need collating at the
    // end of each event
    for (auto & hist: gen_hists) hist.second.collate_event();

    // update event number before checking whether to write
    iev++;
    if (periodic_output_is_due()) standard_output();
  }
  // arrange for final output
  standard_output();
}


//======================================================================
void AnalysisBase::standard_output() {
  cout << "Outputting result after generation of " << iev << " events" << endl;
  
  ofstream ostr(output_filename.c_str());

  ostr << header.str();
  ostr << "# time now = " << cmdline->time_stamp() << endl;
  ostr << "# nev generated = " << iev << endl;

  // normalisation factor to be used throughout
  double norm = weight_factor();

  // write out cross sections
  for(auto label: ordered_labels(xsections)) {
    auto & obj = xsections[label];
    int n = obj.n();
    obj.set_n(iev); // for errors on xsc
    ostr << "# " << label << ": xsc = " << obj.sum()*norm 
         << " +- " << obj.error_on_sum()*norm << " "+_units_string ;
    
    obj.set_n(n); // reset n to be what it was
    ostr << " (n entries = " << obj.n() << ")"
         << endl;
  }

  // write out averages (a little complex because of different scenarios
  // for how we normalise the average)
  for(auto label: ordered_labels(averages)) {
    auto & obj = averages[label];
    int n = obj.n();
    double avnorm;
    const AverageAndError * ref_xsection = 0;
    if (obj.ref_xsection == "" && ! obj.internal_ref) {
      // GPS WARNING: this is only good if we have a total weight
      avnorm = 1.0/total_weight;
      obj.set_n(iev); // for errors on xsc
    } else {
      if (obj.internal_ref) ref_xsection = & obj.ref_xsection_value;
      else                  ref_xsection = & xsections[obj.ref_xsection];
      avnorm = 1.0/ref_xsection->sum();
      obj.set_n(ref_xsection->n()); // for errors
    }
    ostr << "# <" << label << "> = " << obj.sum() * avnorm
         << " +- " << obj.error_on_sum() * avnorm;
    
    obj.set_n(n); // reset n to be what it was
    ostr << " (n entries = " << obj.n();
    if (!ref_xsection) ostr << ", over all events)";
    else if(obj.internal_ref) 
      ostr << ", xsection = " << obj.ref_xsection_value.sum()*norm << " nb)";
    else               ostr << ", wrt "<< obj.ref_xsection << ")";
    ostr << endl;
  }

  // write out normal histograms
  for (const auto & label: ordered_labels(hists)) {
    const auto & obj = hists[label];
    ostr << "# plain_hist:" << label << " [binlo binmid binhi d/dbinvar]" << endl;
    output(obj, &ostr, norm/obj.binsize());
  }

  // output gen_hists
  for (const auto & label: ordered_labels(gen_hists)) {
    const auto & obj = gen_hists[label];
    ostr << "# gen_hist:" << label << " [binlo binmid binhi d/dbinvar err]" << endl;
    obj.write_norm(ostr, norm);
    ostr << endl << endl;
  }

  // output histograms normalised to 1
  for (const auto & label: ordered_labels(norm_hists)) {
    const auto & obj = norm_hists[label];
    ostr << "# norm_hist:" << label << " [binlo binmid binhi d/dbinvar err]" << endl;

    double hist_weight = obj.total_weight();
    ostr << "# histogram total x-sect = " 
         <<  hist_weight * norm << " " << _units_string 
	 << "( " << obj.n_entries() << " entries)"
	 << endl;
    double this_norm = hist_weight != 0.0 ? 1.0/hist_weight : 0;
    
    output(obj, &ostr, this_norm/obj.binsize());
    ostr << endl << endl;
  }

  // output averaging histograms (sometimes called profile histograms)
  for (const auto & label: ordered_labels(avg_hists)) {
    const auto & obj = avg_hists[label];
    ostr << "# avg_hist:" << label << " [binlo binmid binhi; average within each bin; stddev; error on avg; average of squares] " << endl;
    output_noNaN(obj, &ostr);
    ostr << endl << endl;
  }

  // output correlation histograms
  for (const auto & label: ordered_labels(corr_hists)) {
    const auto & obj = corr_hists[label];
    ostr << "# corr_hist:" << label << endl;
    // true for last arg ensures we get more info
    output_noNaN(obj, &ostr, true); 
    ostr << endl << endl;
  }

  // output 2d histograms
  for (const auto & label: ordered_labels(hists_2d)) {
    const auto & obj = hists_2d[label];
    ostr << "# 2d_hist:" << label << " [xlo xmid xhi ylo ymid yhi dN/dxdy] " << endl;
    output(obj, &ostr, norm/obj.u_binsize()/obj.v_binsize());
    ostr << endl << endl;
  }
  
//   output_collection(ostr, norm_hists, "norm_hist:",
//                     [ ](const DefaultHist & hist) {return 1.0/hist.total_weight()/hist.binsize();});
// 
//   output_collection(ostr, avg_hists, "avg_hist:",
//                     [ ](const DefaultAveragingHist & hist) {return 1.0;});
//
}