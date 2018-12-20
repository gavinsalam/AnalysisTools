#ifndef   __ANALYSISFRAMEBASE_HH__
#define   __ANALYSISFRAMEBASE_HH__


#include "SimpleNTuple.hh"
#include "SimpleHist.hh"
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
      this->declare(xmin, xmax, dx);
      _first_call = false;
    }
    this->add_entry(val, weight);
  }

  /// this one can be used for AveragingHist
  void set_lims_add_entry(double xmin, double xmax, double dx, double val, double contents, double weight) {
    if (_first_call) {
      this->declare(xmin, xmax, dx);
      _first_call = false;
    }
    this->add_entry(val, contents, weight);
  }

  /// this one can be used for CorrelationHist
  void set_lims_add_entry(double xmin, double xmax, double dx, double val, double varA, double varB, double weight) {
    if (_first_call) {
      this->declare(xmin, xmax, dx);
      _first_call = false;
    }
    this->add_entry(val, varA, varB, weight);
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

typedef TemplateDefaultHist<SimpleHist>    DefaultHist;
typedef TemplateDefaultHist<AveragingHist> DefaultAveragingHist;
typedef TemplateDefaultHist<CorrelationHist> DefaultCorrelationHist;

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
  void add_entry(double value, double evwgt) {
    *this              += value*evwgt;
    ref_xsection_value += evwgt;
    internal_ref        = true;
  }
  /// the name of the reference cross section that will be used
  /// to normalise this average
  string ref_xsection;
  AverageAndError ref_xsection_value;
  bool   internal_ref;
};

//------------------------------------------------------------

/// Class that helps organise an analysis. The idea is that within a
/// given framework a number of the protected functions will be set up
/// (e.g. to handle a specific type of event record, etc.).
///
/// The base class can handle a variety of aspects
/// 
/// - maps of various kinds of histograms (which can have default
///   binnings, or binnings set up on the fly), and which get
///   automatically output
///
/// - periodic output of results
class AnalysisBase {
public:
  AnalysisBase(CmdLine & cmdline);
  virtual ~AnalysisBase ();

  /// the user may want to set up some of their own parameters
  virtual void user_startup() {};

  /// things that the user will want to do after all the
  /// options have been processed
  virtual void user_post_startup() {};

  /// let the user fill their histograms, etc.
  virtual void analyse_event() {};

  /// any extra output
  virtual void user_output(std::ostream &) {};

protected:

  CmdLine & cmdline;

  /// for things that the framework sets up ahead of time;
  virtual void pre_startup();
  virtual void post_startup();
  virtual void standard_output();

  /// this has to be set up by whatever derives from this class
  //// 
  virtual void event_loop() = 0;

  std::map<std::string,DefaultHist> hists; //< histograms normalised as a differential cross section on output
  std::map<std::string,DefaultHist> norm_hists; //< histograms normalised to have total weight of 1 on output
  std::map<std::string,DefaultAveragingHist> avg_hists;
  std::map<std::string,DefaultCorrelationHist> corr_hists;
  std::map<std::string,SimpleHist2D> hists_2d;
  std::map<std::string,AverageAndError> xsections;
  std::map<std::string,AverageAndErrorWithRef> averages;
  std::map<std::string,NLOHistGeneric> gen_hists;
  typedef std::map<std::string, NLOHistGeneric>::iterator GenHistIt;

  std::ostringstream header;
  std::string output_filename;

  double _total_weight;

  /// this will be used to output units
  std::string _units_string;

};

#endif // __ANALYSISFRAMEBASE_HH__