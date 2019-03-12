#ifndef   __ANALYSISFRAMEBASE_HH__
#define   __ANALYSISFRAMEBASE_HH__


#include <limits>
#include <unordered_map>

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
  AnalysisBase(CmdLine * cmdline_in);

  void run();
  
  virtual ~AnalysisBase () {}

  /// the user may want to set up some of their own parameters
  virtual void user_startup() {}

  /// things that the user will want to do after all the
  /// options have been processed
  virtual void user_post_startup() {}

  /// this generates the event
  virtual void generate_event() {};

  // this returns the weight of the current event
  // (only valid after a call to generate_event)
  virtual double event_weight() {return event_weight_;}
  
  /// let the user fill their histograms, etc.
  virtual void analyse_event() {}

  /// any extra output
  virtual void user_output(std::ostream &) {}

  /// set the default binning for each kind of histogram
  void set_default_binning(double xmin, double xmax, double dx) {
    DefaultHist::set_defaults(xmin, xmax, dx); 
    DefaultAveragingHist::set_defaults(xmin, xmax, dx); 
    DefaultCorrelationHist::set_defaults(xmin, xmax, dx); 
  }

  /// returns the factor that multiplies the sum of weights in order
  /// to get a cross section (or analogue).
  virtual double weight_factor() const {return 1.0/iev;}

  
protected:

  CmdLine * cmdline;

  /// for things that the framework sets up ahead of time;
  virtual void pre_startup() {};
  virtual void post_startup() {};
  virtual void standard_output();


  template<class T> using Collection = std::unordered_map<std::string, T>;

  // scalar type quantities
  Collection<AverageAndError> xsections;
  Collection<AverageAndErrorWithRef> averages;

  // 1D type quantities
  Collection<DefaultHist> hists; //< histograms normalised as a differential cross section on output
  Collection<DefaultHist> norm_hists; //< histograms normalised to have total weight of 1 on output
  Collection<DefaultAveragingHist> avg_hists;
  Collection<DefaultCorrelationHist> corr_hists;
  Collection<NLOHistGeneric> gen_hists;
  typedef Collection<NLOHistGeneric>::iterator GenHistIt;

  // 2d quantities
  Collection<SimpleHist2D> hists_2d;

  
  std::ostringstream header;
  std::string output_filename;

  double event_weight_ = 1.0;
  double total_weight = 0.0;

  /// this will be used to output units
  std::string _units_string;

  CmdLine * cmdline_;

  unsigned long long int iev = 0, nev = 0, output_interval = 1, iev_last_output=0;

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
