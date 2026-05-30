#ifndef   __ANALYSISFRAMEBASE_HH__
#define   __ANALYSISFRAMEBASE_HH__


#include <limits>
#include <unordered_map>

#include "SimpleNTuple.hh"
#include "SimpleHist.hh"
#include "SimpleHistWithError.hh"
#include "SimpleHist2D.hh"
#include "AveragingHist.hh"
#include "CorrelationHist.hh"
#include "AverageAndError.hh"
#include "NLOHistGeneric.hh"

#include "CmdLine.hh"


//----------------------------------------------------------------------
/// a histogram with default limits and bin size
template<class T>
class TemplateDefaultHist : public T {
private:
  bool _first_call;
  static double _lo, _hi, _bin_size;

public:
  /// create a histogram with the default bounds
  TemplateDefaultHist() : T(_lo, _hi, _bin_size), _first_call(true) {}

  /// a way to set non-default limits  and then add an entry
  void set_lims_add_entry(double xmin, double xmax, double dx, double val, double weight = 1.0) {
    if (_first_call) {
      T::declare(xmin, xmax, dx);
      _first_call = false;
    }
    this->add_entry(val, weight);
  }

  /// this one can be used for AveragingHist
  void set_lims_add_entry(double xmin, double xmax, double dx, double val, double contents, double weight) {
    if (_first_call) {
      T::declare(xmin, xmax, dx);
      _first_call = false;
    }
    this->add_entry(val, contents, weight);
  }

  /// this one can be used for CorrelationHist
  void set_lims_add_entry(double xmin, double xmax, double dx, double val, double varA, double varB, double weight) {
    if (_first_call) {
      T::declare(xmin, xmax, dx);
      _first_call = false;
    }
    this->add_entry(val, varA, varB, weight);
  }

  T & declare_once(const Binning & binning) {
    if (_first_call) {
      T::declare(binning);
      _first_call = false;
    }
    return *this;
  }

  T & declare_once(double minv, double maxv, double bin_size) {
    if (_first_call) {
      T::declare(minv,maxv,bin_size);
      _first_call = false;
    }
    return *this;
  }

  TemplateDefaultHist & declare(double minv, double maxv, double bin_size) {
    T::declare(minv,maxv,bin_size);
    _first_call = false;
    return *this;
  }

  TemplateDefaultHist & declare(const Binning & binning) {
    T::declare(binning);
    _first_call = false;
    return *this;
  }

  static void set_defaults(double xmin, double xmax, double dx) {
    _lo = xmin;
    _hi = xmax;
    _bin_size = dx;
  }
  
};

// these definitions are needed to avoid warnings with Apple LLVM version 8.1.0 (clang-802.0.41)
#ifdef __clang__   
template <typename T> double TemplateDefaultHist<T>::_lo = 0.0;
template <typename T> double TemplateDefaultHist<T>::_hi = 1.0;
template <typename T> double TemplateDefaultHist<T>::_bin_size = 0.1;
#endif

typedef TemplateDefaultHist<SimpleHist>          DefaultHist;
typedef TemplateDefaultHist<SimpleHistWithError> DefaultHistWithError;
typedef TemplateDefaultHist<AveragingHist>       DefaultAveragingHist;
typedef TemplateDefaultHist<CorrelationHist>     DefaultCorrelationHist;

//----------------------------------------------------------------------
/// slight variant of average and error, which has the name of
/// a reference cross section for normalisation.
///
/// Just set the variable ref_xsection to be the string that
/// characterises the cross section that is to be used as a reference.
///
class AverageAndErrorWithRef : public AverageAndError {
public:
  AverageAndErrorWithRef() :  internal_ref(false) {}
  inline void add(double x) {
    AverageAndError::add(x);
    add_for_minmax(x);
  }

  inline void add_for_minmax(double x) {
    if (x < _min) _min = x;
    if (x > _max) _max = x;
  }
  
  void add_entry(double value, double evwgt) {
    AverageAndError::add(value*evwgt);
    add_for_minmax(value);
    ref_xsection_value += evwgt;
    internal_ref        = true;
  }
  AverageAndErrorWithRef & set_ref(std::string ref) {
    ref_xsection = ref;
    return *this;
  }
  /// the name of the reference cross section that will be used
  /// to normalise this average
  string ref_xsection;
  AverageAndError ref_xsection_value;
  bool   internal_ref;

    /// return the the minimum value encountered
  double min() const {return _min;}
  /// return the the maximum value encountered
  double max() const {return _max;}
  
private:
  /// internal store of the bounds
  /// (note the _min starts as numeric_limits<double>::max() because then any thing smaller
  /// than that will immediately register as the minimum)
  double _min = std::numeric_limits<double>::max(), _max=-std::numeric_limits<double>::max();
};

/* ********************************************************************/
/// Abstract base class to help organise event generation and analysis
/// and histogram output within a single framework.
///
/// Basic usage:
///
/// - implement \ref generator_routines_section "generator" functions, notably by deriving generate_event(), 
///   and, optionally, pre_startup(), post_startup(), pre_run(),
///   generator_output()
///
/// - implement \ref analysis_routines_section "analysis" functions, notably by deriving user_analyse_event()
///   and, optionally, user_startup(), user_post_startup(), user_analyse_event(),
///   user_output()
///
/// In some settings one might have a class hierarchy 
/// AnalysisBase -> AnalysisWithGenerator -> AnalysisWithGeneratorAndAnalysis,
/// where the middle class implements the generator functions and the final
/// class implements the analysis functions.
///
/// The framework itself provides
///
/// - some basic command-line options (e.g. -nev, -o for the output file)
///
/// - the main event loop 
///
/// - a \ref header ostringstream, e.g. to register/output details of the run
///
/// - maps of histograms of various kinds (which can have
///   default binnings, or binnings set up on the fly)
///
/// - automatic output of histograms, averages, etc. at periodic intervals
///   and at the end of the event loop
///
/// Assuming AnalysisClass is a class derived from AnalysisBase,
/// typical usage
///
/// @code
///   int main(int argc, char * argv[]) {
///     CmdLine cmdline(argc, argv);
///     AnalysisClass analysis(&cmdline);
///     analysis.run();
///     return 0;
///   }
/// @endcode
///
class AnalysisBase {
public:
  AnalysisBase(CmdLine * cmdline_in, const string & default_output_filename = "" );

  /// Once a class derived fromAnalysisBase is fully up, the run()
  /// function carries out all initialisation calls and then runs the
  /// main event loop + any I/O
  void run();
  
  virtual ~AnalysisBase () {}



  /// @anchor generator_routines_section
  /// @name Generator routines
  ///
  /// Routines that derived classes can/should override to provide the
  /// event-generation side of the framework (as opposed to user
  /// analysis). The only routine that **must** be overridden is
  /// generate_event().
  /** @{ **********************************************************/
  
  /// Adapt this routine for generator startup items that you don't want
  /// to have in the constructor and that need the CmdLine. Quite often
  /// this will be empty, with the constructor doing all the work, but
  /// this member can be useful if there is class hierarchy structure
  /// such that you need to be able to access virtual functions during
  /// startup.
  virtual void pre_startup() {}

  /// This gets called after user_startup(), but before
  /// user_post_startup() and pre_run(). CmdLine completeness gets
  /// checked before this call, so do not add any new CmdLine options
  /// here.
  ///
  /// One use for this is to output settings to the `header` ostringstream. 
  virtual void post_startup() {}

  /// this should generate the event and return true if it has done so
  /// successfully; if it returns false, the AnalysisBase main event
  /// loop will not increase the total number of events generated. It
  /// can also set the event_weight_ protected member if it wants to use
  /// a non-unit weight for an event.
  virtual bool generate_event() = 0;

  /// for things a derivative of the framework might do
  /// after event analysis (keeping in mind that anything
  /// to be done beforehand can go into generate_event()).
  /// Can be used, e.g., for analysis timing.
  virtual void post_analyse_event() {}

  /// This is called before the first event is generated, and might
  /// commonly be used for any generator initialisation
  virtual void pre_run() {}

  /// carries out any output that the generator itself wants to send
  /// to the output stream during the periodic and final I/O. 
  /// E.g. may output a summary of any warnings that arose from
  /// the generator
  virtual void generator_output(std::ostream &) {}
  /** @} **************************************************************/



  /// @anchor analysis_routines_section
  /// @name User/Analysis routines
  ///
  /// Routines that derived classes can/should override to provide the
  /// user analysis side of the framework (as opposed to user
  /// generation). The only routine that **must** be overridden is
  /// user_analyse_event(). 
  ///
  /// Typically, user_analyse_event() will want to access the event itself,
  /// as well as the \ref histogram_etc_access "histogram and cross-section variables".
  /** @{ *************************************************************/
  /// a user use this to set up some of their own parameters
  virtual void user_startup() {}

  /// things that the user will want to do after all the
  /// options have been processed
  virtual void user_post_startup() {}

  /// let the user fill their histograms, etc.
  virtual void user_analyse_event() {}

  /// any extra output
  virtual void user_output(std::ostream &) {}
  /** @} *************************************************************/



  /// @anchor useful_variables
  /// @name Useful variables and member functions for derived classes to access
  /** @{ *************************************************************/

  /// return the weight of the current event (only valid after a
  /// call to generate_event); already has a sensible implementation, but
  /// can be overridden if needed.
  virtual double event_weight() {return event_weight_;}

  /// current event index
  unsigned long long int iev = 0;
  /// requested total number of events
  unsigned long long int nev = 0;
  
  /// set the default binning to some common value the main 1D histograms
  void set_default_binning(double xmin, double xmax, double dx) {
    DefaultHist::set_defaults(xmin, xmax, dx);
    DefaultHistWithError::set_defaults(xmin, xmax, dx);
    DefaultAveragingHist::set_defaults(xmin, xmax, dx); 
    DefaultCorrelationHist::set_defaults(xmin, xmax, dx); 
  }

  /// @brief set the output precision for histograms and averages;
  /// @param prec the number of digits (as in ostream::precision)
  void set_output_precision(int prec) {
    output_precision_ = prec;
  }  

  /// @brief  an ostringstream object for the header
  /// @details This is to be used by generator/user routines to output header information,
  /// while the current state is accessed with header.str()
  std::ostringstream header;
  /// @brief  the output filename, e..g as set by the -o command-line option
  std::string output_filename;

  /** @} *************************************************************/


  /// @anchor histogram_etc_access
  /// @name Access to histograms, averages, etc.
  ///
  /// A range of "Collection" (`unordered_map<std::string, T>`) objects
  /// for storing cross sections, histograms, etc. These are
  /// automatically written to the output with the correct normalisation
  /// at periodic intervals and at the end of the event loop.
  ///
  /// For example
  /// @code
  ///   xsections["my_cross_section"] += event_weight();
  ///   hists["my_histogram"].add_entry(value, event_weight());
  /// @endcode
  /** @{ >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

  template<class T> using Collection = std::unordered_map<std::string, T>;

  // scalar type quantities
  Collection<AverageAndError> xsections;
  Collection<AverageAndErrorWithRef> averages;

  // 1D type quantities
  /// histograms normalised as a differential cross section on output
  Collection<DefaultHist> hists;
  /// similar to hists, but one outputs the cumulative histogram as well
  Collection<DefaultHist> cumul_hists;
  /// histograms normalised as a differential cross section on output
  Collection<DefaultHistWithError> hists_err;
  /// unit normalised histograms with errors
  Collection<DefaultHistWithError> norm_hists_err;
  /// similar to hists, but one outputs the cumulative histogram as well
  Collection<DefaultHistWithError> cumul_hists_err;
  /// histograms normalised to have total weight of 1 on output
  Collection<DefaultHist> norm_hists; 
  Collection<DefaultAveragingHist> avg_hists;
  Collection<DefaultCorrelationHist> corr_hists;
  Collection<NLOHistGeneric> gen_hists;
  typedef Collection<NLOHistGeneric>::iterator GenHistIt;

  // 2d quantities
  /// 2d-histograms, where all bin edges get written out
  Collection<SimpleHist2D> hists_2d;
  /// 2d-histograms, where only bin midpoint coordinates get written out
  Collection<SimpleHist2D> hists_2d_compact;

  /** @} *************************************************************/


  /// @anchor external_event_loop_section
  /// @name Members for use with external event loops
  ///
  /// Functions for users who wish to use a class derived from
  /// AnalysisBase in a more manual way, handling initialisation and the
  /// event loop themselves.
  /** @{ **********************************************************/
  /// carry out all initialisation calls
  void do_all_init();

  /// generate the next event and return true if the event has been
  /// successfully generated. Note that \ref iev index is incremented
  /// even if the event is not successful and the total number of events attempted/generated is iev+1.
  ///
  /// Note: the requested number of events, \ref nev, is ignored, i.e.
  /// iev+1 can exceed nev. Output is not triggered by this function.
  ///
  /// Relative to calling generate_event(), the main differences are that 
  /// this routine makes sure initialisation has been carried out and that
  /// the iev index is up to date.
  bool next_event() {
    if (!_init_done) do_all_init();
    if (_first_event_done) ++iev;
    _first_event_done = true;
    bool output =  generate_event();
    total_weight += event_weight();
    return output;
  }
  /** @} **********************************************************/



  /// @anchor internal_variables
  /// @name Variables and members mainly intended for internal use
  /** @{ *************************************************************/

  /// return the factor that multiplies the sum of weights in order
  /// to get a cross section (or analogue).
  virtual double weight_factor() const {return 1.0/iev;}
  /// return the dimensions of the cross section (or analogue)
  virtual std::string units_string() const {return "arb.units";}

  /// return true if it's time for periodic output (established based on
  /// the event number); this will not check if periodic output actually
  /// occurred, but still updates internal counters as if it had
  bool periodic_output_is_due() {
    if (iev - iev_last_output < output_interval) return false;
    iev_last_output = iev;
    if (output_interval * (1.0/iev) < 0.05000000001) {
      if (time_since_last_write() < _max_time_between_writes) output_interval *= 2;
    }
    return (time_since_last_write() >= _min_time_between_writes);
    //return true;
  }

  /// return the time (in whole s) since the last write 
  double time_since_last_write() const {
    return cmdline->time_elapsed_since_start() - _time_elapsed_at_last_write;
  }
  /** @} *************************************************************/
  

protected:

  CmdLine * cmdline;


  virtual void standard_output();
  
  /// In some generation schemes, one may make many attempts at phase
  /// space generation for a single successful event; to take this into
  /// account in error estimates, we reset the n in xsection and average
  /// and error objects to the number returned by this function. 
  ///
  /// By default the function just returns iev, as is good, e.g. if 
  /// the normalising cross section is known exactly. 
  /// 
  /// However derived classes may override this function. For example if
  /// the generator has a fixed overhead for a Born cross section and
  /// has made 1000 trials in order to successfully generate 100 events
  /// then by returning 1000 rather than 100, this function enables
  /// output of cross sections to correctly generate the relevant errors.
  ///
  /// Note this is only relevant for the error; the actual normalisation
  /// should be dealt with by overriding the weight_factor() function 
  /// above.
  virtual uint64_t effective_iev_attempts() const {return iev;}

  
  /// current even weight (only valid after a call to generate_event)
  double event_weight_ = 1.0;
  /// sum of weights of all successfully generated events
  double total_weight = 0.0;

  /// the current interval (in number of events) between writes to output
  unsigned long long int output_interval = 10;
  /// the iev index, plus one, at the time periodic_output_is_due()
  /// last returned true
  unsigned long long int iev_last_output=0;
  /// the maximum run-time in seconds
  double max_time_s = -1.0;
  /// the elapsed time at which the last write occurred (in whole s)
  double _time_elapsed_at_last_write = 0.0;
  /// the minimum time between writes (in s), only whole seconds registered
  double _min_time_between_writes = 1.0;
  /// the approx maximum time between writes (in s; write frequency stops doubling if interval exceeds this)
  double _max_time_between_writes = 900;

  template<class T> vector<string> ordered_labels(Collection<T> & collection) {
    std::vector<string> labels;
    // first get an alphabetically ordered list of labels, regardless
    // of whether the collection has an underlying alphabetical order (map)
    // or not (unordered_map)
    for(const auto & item: collection) {
      labels.push_back(item.first);
    }
    std::sort(labels.begin(), labels.end());
    return labels;
  }
  
  int output_precision_ = - 1;

  bool _init_done = false;
  bool _first_event_done = false;

//   /// output a given collection, in alphabetical order of the labels
//   template<class T, class U> std::ostream & output_collection(std::ostream & ostr,
//                                                               Collection<T> & collection,
//                                                               const std::string & prefix,
//                                                               U norm_function) {
//     std::vector<string> labels;
//     // first get an alphabetically ordered list of labels, regardless
//     // of whether the collection has an underlying alphabetical order (map)
//     // or not (unordered_map)
//     for(const auto & item: collection) {
//       labels.push_back(item.first);
//     }
//     std::sort(labels.begin(), labels.end());
//     for (const auto & label: labels) {
//       ostr << "# " << prefix << ":" << label << endl;
//       const auto & item = collection[label];
//       auto norm = norm_function(item);
//       output(item, &ostr, norm);
//     }
//     return ostr;
//   };
};

#endif // __ANALYSISFRAMEBASE_HH__
