#ifndef __CORRELATIONHIST_HH__
#define __CORRELATIONHIST_HH__

#include "SimpleHist.hh"

/// class that provides an "histogram" of correlation coefficients
class CorrelationHist {
public:

  CorrelationHist() : _first_call(true){}

  CorrelationHist(double minv, double maxv, int n) {
    declare(minv, maxv, n);
  }

  CorrelationHist(double minv, double maxv, unsigned int n) {
    declare(minv, maxv, n);
  }

  CorrelationHist(double minv, double maxv, double bin_size) {
    declare(minv, maxv, int(0.5+(maxv-minv)/bin_size));
  }

  // declare (or redeclare) the histogram
  void declare(double minv, double maxv, double bin_size) {
    declare(minv, maxv, int(0.5+(maxv-minv)/bin_size));
  }

  void declare(double minv, double maxv, int n) {
    declare(minv, maxv, unsigned(n));
  }

  // declare (or redeclare) the histogram
  void declare(double minv, double maxv, unsigned int n) {
    _sumx.declare(minv, maxv, n);
    _sumy.declare(minv, maxv, n);
    _sumxx.declare(minv, maxv, n);
    _sumxy.declare(minv, maxv, n);
    _sumyy.declare(minv, maxv, n);
    _weights.declare(minv, maxv, n);
    _nentries.declare(minv, maxv, n);
    _first_call = false;
  }

  void add_entry(double v, double x, double y, double weight = 1.0) {
    _sumx.add_entry (v, x * weight);
    _sumy.add_entry (v, y * weight);
    _sumxx.add_entry(v, x*x * weight);
    _sumxy.add_entry(v, x*y * weight);
    _sumyy.add_entry(v, y*y * weight);
    _weights.add_entry(v, weight);
    _nentries.add_entry(v, 1.0);
  };

  void set_lims_add_entry(double vmin, double vmax, double dv, double v, double x, double y, double weight) {
    if (_first_call) {
      declare(vmin, vmax, dv);
      _first_call = false;
    }
    add_entry(v, x, y, weight);
  }

  double binlo (int i) const {return _sumx.binlo(i);}
  double binhi (int i) const {return _sumx.binhi(i);}
  double binmid(int i) const {return _sumx.binmid(i);}
  double binsize()     const {return _sumx.binsize();}
  unsigned int bin(double v) const {return _sumx.bin(v);}

  double min() const {return _sumx.min();};
  double max() const {return _sumx.max();};
  /// returns the size of the histogram proper
  unsigned int size() const {return _sumx.size();};
  /// returns the size of the histogram plus outflow bin
  unsigned int outflow_size() const {return _sumx.outflow_size();};
  
  /// return the average of the current bin
  double correlation(int i) const {
    if (_nentries[i]<=1) return 0.0;// we need at the very least 2 points
    double avgx  = _sumx[i] /_weights[i];
    double avgy  = _sumy[i] /_weights[i];
    double avgxx = _sumxx[i]/_weights[i];
    double avgxy = _sumxy[i]/_weights[i];
    double avgyy = _sumyy[i]/_weights[i];
    // the line below should be multiplied by (n-1)/(n-2) to get an
    // unbiased estimate of the correlation when using a finite "sample"
    return (avgxy-avgx*avgy)/sqrt(std::abs(avgxx-avgx*avgx)*std::abs(avgyy-avgy*avgy));
  }

  double avgx (int i) const {return _sumx[i]/_weights[i];}
  double avgy (int i) const {return _sumy[i]/_weights[i];}
  double avgxy(int i) const {return _sumxy[i]/_weights[i];}
  double varx (int i) const {return _sumxx[i]/_weights[i] - pow(avgx(i),2);}
  double vary (int i) const {return _sumyy[i]/_weights[i] - pow(avgy(i),2);}

private:
  SimpleHist _sumx, _sumy, _sumxx, _sumxy, _sumyy, _weights, _nentries;
  bool _first_call;
};

//----------------------------------------------------------------------
inline void output_noNaN(const CorrelationHist & hist, 
                         std::ostream * ostr = (&std::cout),
                         bool details = false) {
  
  *ostr << "# binlo binmid binhi correl";
  if (details) *ostr << " avgx avgy avgxy varx vary ";
  *ostr << std::endl;
  for (unsigned i = 0; i < hist.size(); i++) {
    *ostr << hist.binlo(i)  << " " 
          << hist.binmid(i) << " "
          << hist.binhi(i) << " "
          << hist.correlation(i);
    if (details) {
      *ostr 
        << " " << hist.avgx (i)
        << " " << hist.avgy (i)
        << " " << hist.avgxy(i)
        << " " << hist.varx (i)
        << " " << hist.vary (i);
    }
    *ostr << std::endl;
  }
}


#endif //  __AVERAGINGHIST_HH__
