#ifndef __BINNING_HH__
#define __BINNING_HH__
#include<cassert>

class Binning {
public:
  Binning()  : _declare_calls_init(true) {}

  Binning(double minv, double maxv, unsigned int n) {
    declare(minv, maxv, n);
    /// so far this was false; but by setting true, then calls
    /// to declare(...) from derived classes will call that
    /// derived class's _init() function
    _declare_calls_init = true;
  }

  // NB: the int case is needed because if one passes an
  // int then compiler doesn't know whether to convert
  // to unsigned or double.
  Binning(double minv, double maxv, int n) {
    declare(minv, maxv, unsigned(n));
    _declare_calls_init = true;
  }

  Binning(double minv, double maxv, double bin_size) {
    declare(minv, maxv, bin_size);
    _declare_calls_init = true;
  }

  // declare (or redeclare) the histogram
  void declare(double minv, double maxv, double bin_size) {
    declare(minv, maxv, int(0.5+std::abs((maxv-minv)/bin_size)));
  }

  void declare(double minv, double maxv, int n) {
    declare(minv, maxv, unsigned(n));
  }

  // declare (or redeclare) the histogram
  void declare(double minv, double maxv, unsigned int n) {
    _minv = minv; _maxv = maxv; _dv = (maxv-minv)/n;
    _size = n;
    _outflow_size = n+2;
    _init();
  }

  void declare(const Binning & binning) {
    *this = binning;
    _init();
  }

  // this does nothing in this class, but is triggered when binning is
  // declared and allows derived classes to do something (without
  // their having to reimplement all the declare routines)
  virtual void _init() {}
  
  /// returns the lower extent of the binning
  double min() const {return _minv;}
  /// returns the upper extent of the binning
  double max() const {return _maxv;}

  /// returns the size of the histogram proper (excluding outflow bins)
  unsigned int size() const {assert (_outflow_size > 1); return _size;}
  /// returns the size of the histogram plus outflow bin
  unsigned int outflow_size() const {assert (_outflow_size > 1); return _outflow_size;}

  /// returns the lower edge of bin i
  double binlo (int i) const {return i*_dv + _minv;};
  /// returns the upper edge of bin i
  double binhi (int i) const {
    // special condition for bin i == _size-1 avoids rounding errors on the upper edge
    // (in PanScales these were causing problems in validation runs for the last bin of the histogram which ended at 0)
    return i+1 == int(_size) ? _maxv : (i+1)*_dv + _minv; };
  /// returns the arithmetic midpoint of bin i
  double binmid(int i) const {return (i+0.5)*_dv + _minv;};
  /// returns the binsize
  double binsize()     const {return _dv;};

  /// returns the bin index for a value v
  unsigned int bin(double v) const {
    // divide by _dv before checking under/overflow, because 
    // the sign of _dv (which is not prescribed) affects whether
    // it is under or overflow. BUT: this does can problems if 
    // trapping floating point overflow and we have v that
    // is just at the edge of the numerical domain
    double v_minus_minv_over_dv = (v-_minv)/_dv;

    // handle underflow case
    if (v_minus_minv_over_dv < 0.0) return underflow_bin();
    // handle overflow case
    if (v_minus_minv_over_dv >= size()) return overflow_bin();

    // otherwise just return a bin index, which is bound to be in
    // range given the above two tests
    return unsigned(v_minus_minv_over_dv);
  }
  
  /// returns the index of the underflow bin
  unsigned int underflow_bin() const {return size();}

  /// returns the index of the overflow bin
  unsigned int overflow_bin() const {return size()+1;}

  /// virtual dummy destructor
  virtual ~Binning() {}
  
protected:
  double _minv, _maxv, _dv;
  unsigned _size=0, _outflow_size=0;

  
  /// while the Binning is being constructed, one does not want
  /// declare(...) to calls _init (which is overloaded virtual in
  /// derived classes); once the Binning class is constructed, then
  /// this changes so that calls to declare() do call _init (including
  /// up the chain in the derived class).
  bool _declare_calls_init = false;
};

#endif // __BINNING_HH__
