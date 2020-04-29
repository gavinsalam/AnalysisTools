#ifndef __SIMPLEHISTWITHERROR_HH__
#include "SimpleHist.hh"
class SimpleHistWithError : public SimpleHist {
public:
  SimpleHistWithError() {}
  SimpleHistWithError(const Binning & binning) : SimpleHist(binning) {_init();}
  SimpleHistWithError(double minv, double maxv, unsigned int n) : SimpleHist(minv,maxv,n) {_init();}

  // NB: the int case is needed because if one passes an
  // int then compiler doesn't know whether to convert
  // to unsigned or double.
  SimpleHistWithError(double minv, double maxv, int n) : SimpleHist(minv,maxv,n) {
    _init();
  }

  SimpleHistWithError(double minv, double maxv, double bin_size) : SimpleHist(minv,maxv,bin_size)  {
    _init();
  }

  void _init() override {
    SimpleHist::_init();
    _weights_sumsqr.resize(outflow_size());
  }
  void reset() {
    SimpleHist::reset();
    _weights_sumsqr = 0.0;
  }
  /// returns the sum of squared weights in the bin
  double & sumsqr(unsigned i)  {return _weights_sumsqr[i];}

  /// returns the sum of squared weights in the bin
  const double & sumsqr(unsigned i) const {return _weights_sumsqr[i];}

  /// returns the error on the bin's contents
  double error(unsigned i) const {return error_calc(_weights[i], _weights_sumsqr[i]);}

  // Operations with constants ---------------------------------------
  SimpleHistWithError & operator*=(double fact) {
    SimpleHist::operator*=(fact);
    double factsqr = fact*fact;
    for (unsigned i = 0; i < outflow_size(); i++) sumsqr(i) *= factsqr;
    return *this;
  }
  SimpleHistWithError & operator/=(double fact) {
    *this *= 1.0/fact;
    return *this;
  }

  
  //--------- IO operations -------------------
  /// output the histogram, including the header, as a differential
  /// histogram d/dv. Outflow bins don't get any bin width normalisaion
  /// The prefix is used for lines not intended to be used in normal plots
  std::ostream & output(std::ostream & ostr, const std::string & prefix = "# ") {
    output_total_and_mean(ostr, prefix);
    ostr << prefix << "cols: vlo vmid vhi hist[v] err[v]" << std::endl;
    ostr << prefix << "under flow bin " << underflow() << " +- " << error(underflow_bin()) << std::endl;
    for (unsigned i = 0; i < size() ; i++) {
      ostr << binlo(i)  << " " 
           << binmid(i) << " "
           << binhi(i) << " "
           << (*this)[i] << " " 
           << error(i)
           << std::endl;
    }
    ostr << prefix << "over flow bin " << overflow() << " +- " << error(overflow_bin()) << std::endl;
    return ostr;
  }

  std::ostream & output_diff(std::ostream & ostr, const std::string & prefix = "# ") {
    output_total_and_mean(ostr, prefix);
    ostr << prefix << "cols: vlo vmid vhi dhist/dv derr/dv (outflow not normalised)" << std::endl;
    ostr << prefix << "under flow bin " << underflow() << " +- " << error(underflow_bin()) << std::endl;
    for (unsigned i = 0; i < size() ; i++) {
      ostr << binlo(i)  << " " 
           << binmid(i) << " "
           << binhi(i) << " "
           << (*this)[i] / (binhi(i) - binlo(i)) << " " 
           << error(i) / (binhi(i) - binlo(i))
           << std::endl;
    }
    ostr << prefix << "over flow bin " << overflow() << " +- " << error(overflow_bin()) << std::endl;
    return ostr;
  }
  
  /// output the cumulative histogram.
  /// The prefix is used for lines not intended to be used in normal plots
  std::ostream & output_cumul(std::ostream & ostr, const std::string & prefix = "# ") {
    output_total_and_mean(ostr, prefix);
    double cumul   = underflow();
    double cumulsq = sumsqr(underflow_bin());
    ostr << prefix << "cols: v hist_integral_up_to_v err" << std::endl;
    ostr << min() << " " << cumul << " " << error_calc(cumul, cumulsq) << std::endl;
    for (unsigned i = 0; i < size(); i++) {
      cumul   += (*this)[i];
      cumulsq += sumsqr(i);
      ostr << binhi(i) << " " << cumul << " " << error_calc(cumul, cumulsq) << std::endl;
    }
    cumul   += overflow();
    cumulsq += sumsqr(overflow_bin());
    ostr << prefix << "with_overflow " << cumul << " +- " << error_calc(cumul, cumulsq) << std::endl;
    return ostr;
  }
  
  friend SimpleHistWithError operator*(const SimpleHistWithError & hist, double fact);
  friend SimpleHistWithError operator/(const SimpleHistWithError & hist, double fact);


protected:
  void _add_entry_ibin(double v, unsigned ibin, double weight) override {
    SimpleHist::_add_entry_ibin(v, ibin, weight);
    _weights_sumsqr[ibin] += weight*weight;
  }

  double error_calc(double sum, double sumsq) const {
    return std::sqrt(std::abs(sumsq - sum*sum/n_entries()));
  }
  std::valarray<double> _weights_sumsqr;
};

inline SimpleHistWithError operator*(double fact, const SimpleHistWithError & hist) {
  return hist*fact;
}
inline SimpleHistWithError operator*(const SimpleHistWithError & hist, double fact) {
  SimpleHistWithError result(hist);
  result *= fact;
  return result;
}
inline SimpleHistWithError operator/(const SimpleHistWithError & hist, double fact) {
  SimpleHistWithError result(hist);
  result *= 1.0/fact;
  return result;
}

#define __SIMPLEHISTWITHERROR_HH__
#endif // __SIMPLEHISTWITHERROR_HH__
