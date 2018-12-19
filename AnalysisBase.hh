#ifndef   __ANALYSISFRAMEBASE_HH__
#define   __ANALYSISFRAMEBASE_HH__

#include "SimpleNTuple.hh"
#include "SimpleHist.hh"
#include "SimpleHist2D.hh"
#include "AveragingHist.hh"
#include "CorrelationHist.hh"
#include "AverageAndError.hh"
#include "NLOHistGeneric.hh"

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



#endif // __ANALYSISFRAMEBASE_HH__