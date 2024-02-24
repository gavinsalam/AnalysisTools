#include "AnalysisBase.hh"
#include <fstream>
#include <exception>

using namespace std;

template<> double DefaultHist::_lo =   0.0;
template<> double DefaultHist::_hi = 100.0;
template<> double DefaultHist::_bin_size = 1.0;

template<> double DefaultHistWithError::_lo =   0.0;
template<> double DefaultHistWithError::_hi = 100.0;
template<> double DefaultHistWithError::_bin_size = 1.0;

template<> double DefaultAveragingHist::_lo =   0.0;
template<> double DefaultAveragingHist::_hi = 100.0;
template<> double DefaultAveragingHist::_bin_size = 1.0;

template<> double DefaultCorrelationHist::_lo =   0.0;
template<> double DefaultCorrelationHist::_hi = 100.0;
template<> double DefaultCorrelationHist::_bin_size = 1.0;


//======================================================================
AnalysisBase::AnalysisBase(CmdLine * cmdline_in,
                           const string & default_output_filename) :
                           cmdline(cmdline_in) {
  cmdline->start_section("Run steering");
  
  nev = cmdline->value<double>("-nev", 1e2).help("number of events to generate").argname("nev");
  max_time_s = cmdline->value<double>("-max-time",-1.0).help("Maximum time (in seconds) that the code should run for");
  output_interval = cmdline->value<double>("-output-interval", output_interval).
                      help("sets the initial output interval (which is then progressively increased)");

  if(default_output_filename.empty()) {
    output_filename = cmdline->any_value<string>({"-o","-out"})
                              .help("filename for output").argname("filename");
  } else {
    output_filename = cmdline->any_value<string>({"-o","-out"}, default_output_filename)
                              .help("filename for output").argname("filename");
  }
  
  bool user_git_info = cmdline->value_bool("-git-info", true)
         .help("include runtime git status in the header() output. "
               "Use '-no-git-info' or '-git-info no' to disable, e.g. on systems where git repo access is slow").argname("yes/no").no_dump();
  cmdline->set_git_info_enabled(user_git_info);

  cmdline->end_section("Run steering");

  header << cmdline->header();

}

//======================================================================
void AnalysisBase::run() {
  pre_startup();
  user_startup();
  
  cmdline->start_section("Run steering");
  auto dump_opt = cmdline->optional_value<string>("-dump-argfile")
              .help("dump the command line arguments to a file (which can be reread with the -argfile option)")
              .argname("filename").no_dump();
  cmdline->end_section("Run steering");

  cmdline->assert_all_options_used();

  if (dump_opt.present()) {
    ofstream ostr(dump_opt());
    if (!ostr.good()) throw runtime_error("Could not write to "+dump_opt());
    ostr << cmdline->dump() << endl;
  }

  post_startup();
  user_post_startup();

  pre_run();
  header << "# time after pre_run initialisation = " << cmdline->time_stamp() ;
  double time_elapsed = cmdline->time_elapsed_since_start();
  header << ", elapsed = " << time_elapsed << " s = " << time_elapsed/60.0 << "m = " << time_elapsed/3600.0 << " h";
  header << endl;
  
  /// should this go into a separate event_loop routine?
  for (iev = 0; iev < nev; ) {

    bool success = generate_event();
    // what do we do on failure: bail out, or continue round the loop?
    // As of 2020-01-23 (GPS), bail out for this event (other than
    // potential periodic output and iev update), but continue
    // with subsequent events
    if (success) {
      total_weight += event_weight();
  
      user_analyse_event();
  
      // certain histogram types need collating at the
      // end of each event
      for (auto & hist: gen_hists) hist.second.collate_event();
  
      // do anything else that's needed at the end of the analysis
      post_analyse_event();
    }

    // update event number before checking whether to write
    iev++;
    if (periodic_output_is_due()) standard_output();

    // allow an escape route if we've exceeded the allowed time
    if (max_time_s > 0.0 && cmdline->time_elapsed_since_start() > max_time_s) break;
  }
  // arrange for final output
  standard_output();
}


//======================================================================
void AnalysisBase::standard_output() {
  cout << "Outputting result after generation of " << iev << " events" << endl;
  _time_elapsed_at_last_write = cmdline->time_elapsed_since_start();
  
  ofstream ostr(output_filename.c_str());
  if (!ostr.good()) throw runtime_error("Could not write to "+output_filename);
  if (output_precision_ > 0) ostr.precision(output_precision_);

  ostr << header.str();
  ostr << "#---------------------------------------------------- " << endl;
  ostr << "# time now = " << cmdline->time_stamp();
  std::streamsize prec = ostr.precision(3);
  double time_elapsed = cmdline->time_elapsed_since_start();
  ostr << ", elapsed = " << time_elapsed << " s = " << time_elapsed/60.0 <<" m = " << time_elapsed/3600.0 << " h";
  ostr << endl;
  ostr.precision(prec);
  ostr << "# nev = " << iev << endl;
  ostr << "#---------------------------------------------------- " << endl;

  // normalisation factor to be used throughout
  double norm = weight_factor();

  // write out cross sections
  for(auto label: ordered_labels(xsections)) {
    auto & obj = xsections[label];
    int n = obj.n();
    obj.set_n(iev); // for errors on xsc
    ostr << "# " << label << ": xsc = " << obj.sum()*norm 
         << " +- " << obj.error_on_sum()*norm << " "+units_string() ;
    
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
      ostr << ", xsection = " << obj.ref_xsection_value.sum()*norm << " " << units_string() << ")";
    else               ostr << ", wrt "<< obj.ref_xsection << ")";
    ostr << ", range=[" << obj.min() << " -- " << obj.max() << "]";
    ostr << endl;
  }

  ostr << "#---------------------------------------------------- " << endl;
  
  // write out normal histograms
  for (const auto & label: ordered_labels(hists)) {
    const auto & obj = hists[label];
    ostr << "# diff_hist:" << label << endl;
    (norm*obj).output_diff(ostr) << endl << endl;
  }

  for (const auto & label: ordered_labels(cumul_hists)) {
    const auto & obj = cumul_hists[label];
    ostr << "# diff_hist:" << label << endl;
    (norm*obj).output_diff(ostr) << endl << endl;

    ostr << "# cumul_hist:" << label << endl;
    (norm*obj).output_cumul(ostr) << endl << endl;
  }

  // write out normal histograms with errors and cumulative variants
  for (const auto & label: ordered_labels(hists_err)) {
    const auto & obj = hists_err[label];
    ostr << "# diff_hist:" << label << endl;
    ostr << "# ecol = 5" << endl;
    (norm*obj).output_diff(ostr) << endl << endl;
  }

  for (const auto & label: ordered_labels(cumul_hists_err)) {
    const auto & obj = cumul_hists_err[label];
    ostr << "# diff_hist:" << label << endl;
    ostr << "# ecol = 5" << endl;
    (norm*obj).output_diff(ostr) << endl << endl;

    ostr << "# cumul_hist:" << label << endl;
    ostr << "# ecol = 3" << endl;
    (norm*obj).output_cumul(ostr) << endl << endl;
  }

  // output gen_hists
  for (const auto & label: ordered_labels(gen_hists)) {
    const auto & obj = gen_hists[label];
    ostr << "# diff_gen_hist:" << label << " [binlo binmid binhi d/dbinvar err]" << endl;
    obj.write_norm(ostr, norm);
    ostr << endl << endl;
  }

  // output histograms normalised to 1
  for (const auto & label: ordered_labels(norm_hists_err)) {
    const auto & obj = norm_hists_err[label];
    ostr << "# norm_diff_hist:" << label << " [binlo binmid binhi d/dbinvar err]" << endl;
    ostr << "# ecol = 5" << endl;

    double hist_weight = obj.total_weight();
    ostr << "# histogram total x-sect = " 
         <<  hist_weight * norm << " " << units_string() 
         << "( " << obj.n_entries() << " entries)"
         << endl;
    double this_norm = hist_weight != 0.0 ? 1.0/hist_weight : 0;
    (this_norm*obj).output_diff(ostr) << endl << endl;
  }

  for (const auto & label: ordered_labels(norm_hists)) {
    const auto & obj = norm_hists[label];
    ostr << "# norm_hist:" << label << " [binlo binmid binhi d/dbinvar err]" << endl;

    double hist_weight = obj.total_weight();
    ostr << "# histogram total x-sect = " 
         <<  hist_weight * norm << " " << units_string() 
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
    ostr << "# ecol = 6" << endl;
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

  // output normal 2d histograms
  for (const auto & label: ordered_labels(hists_2d)) {
    const auto & obj = hists_2d[label];
    ostr << "# 2d_hist:" << label << " [xlo xmid xhi ylo ymid yhi dN/dxdy] " << endl;
    output(obj, &ostr, norm/obj.u_binsize()/obj.v_binsize());
    ostr << endl << endl;
  }

  // output compact 2d histograms
  for (const auto & label: ordered_labels(hists_2d_compact)) {
    const auto & obj = hists_2d_compact[label];
    ostr << "# 2d_hist_compact:" << label << " [xmid ymid dN/dxdy] " << endl;
    output_compact(obj, &ostr, norm/obj.u_binsize()/obj.v_binsize());
    ostr << endl << endl;
  }

  // provide any additional output, for example warnings (maybe
  // warnings should be incorporated into AnalysisTools?)
  user_output(ostr);
  generator_output(ostr);
}
