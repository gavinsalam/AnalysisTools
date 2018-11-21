#ifndef __SIMPLEHIST2D_HH__
#define __SIMPLEHIST2D_HH__

#include<valarray>
#include<string>
#include<cmath>
#include<iostream>
#include<cassert>

class SimpleHist2D {
public:
  SimpleHist2D() {};
  SimpleHist2D(double minu, double maxu, int nu,
	       double minv, double maxv, int nv) {
    declare(minu, maxu, unsigned(nu), minv, maxv, unsigned(nv));
  }

  SimpleHist2D(double minu, double maxu, unsigned int nu,
	       double minv, double maxv, unsigned int nv) {
    declare(minu, maxu, nu, minv, maxv, nv);
  }

  SimpleHist2D(double minu, double maxu, double bin_size_u,
	       double minv, double maxv, double bin_size_v) {
    declare(minu, maxu, bin_size_u, minv, maxv, bin_size_v);
  }

  // declare (or redeclare) the histogram
  void declare(double minu, double maxu, double bin_size_u,
	       double minv, double maxv, double bin_size_v) {
    declare(minu, maxu, int(0.5+(maxu-minu)/bin_size_u),
	    minv, maxv, int(0.5+(maxv-minv)/bin_size_v));
  }

  //void declare(double minv, double maxv, int n) {
  //  declare(minv, maxv, unsigned(n));
  //}

  void declare(double minu, double maxu, int nu,
	       double minv, double maxv, int nv) {
    declare(minu, maxu, unsigned(nu), minv, maxv, unsigned(nv));
  }

  // declare (or redeclare) the histogram
  void declare(double minu, double maxu, unsigned int nu,
	       double minv, double maxv, unsigned int nv) {
    _minu = minu; _maxu = maxu; _du = (maxu-minu)/nu; _nu = nu; 
    _minv = minv; _maxv = maxv; _dv = (maxv-minv)/nv; _nv = nv; 
    _weights.resize(nu*nv+1);
    _weights = 0.0;
    //_weight_v = 0.0;
    //_weight_vsq = 0.0;
    _total_weight = 0.0;
    _have_total = false;
  }

  double u_min() const {return _minu;};
  double u_max() const {return _maxu;};
  unsigned nu () const {return _nu;};
  double v_min() const {return _minv;};
  double v_max() const {return _maxv;};
  unsigned nv () const {return _nv;};
  /// returns the size of the histogram proper
  unsigned int size() const {return _weights.size()-1;};
  /// returns the size of the histogram plus outflow bin
  unsigned int outflow_size() const {return _weights.size();};


  /// return the actual bin for storing iu,iv, alternatively the
  /// outflow bin
  unsigned int getbin(int iu, int iv) const {
    if (iu >= 0 && iu < int(_nu) && iv >= 0 && iv < int(_nv)) {
      return iu*_nv + iv;
    } else {
      return size();
    }
  }

    double & operator()(int iu, int iv) {_have_total = false; 
    return _weights[getbin(iu,iv)];}
    const double & operator()(int iu, int iv) const {
    return _weights[getbin(iu,iv)];}

  double & operator[](int i) {_have_total = false; return _weights[i];}
  const double & operator[](int i) const {return _weights[i];}
  
  /// returns the outflow bin
  double & outflow() {return _weights[size()];};
  const double & outflow() const {return _weights[size()];};

  double u_binlo (int i) const {return i*_du + _minu;};
  double u_binhi (int i) const {return (i+1)*_du + _minu;};
  double u_binmid(int i) const {return (i+0.5)*_du + _minu;};
  double u_binsize()     const {return _du;};

  double v_binlo (int i) const {return i*_dv + _minv;};
  double v_binhi (int i) const {return (i+1)*_dv + _minv;};
  double v_binmid(int i) const {return (i+0.5)*_dv + _minv;};
  double v_binsize()     const {return _dv;};

  unsigned int u_bin(double u) const {
    if (u >= _minu && u < _maxu) {
      int i = int((u-_minu)/_du); 
      if (i >= 0 && i < int(_nu)) {return unsigned(i);} 
    }
    // otherwise...
    return _nu;
  }
  unsigned int v_bin(double v) const {
    if (v >= _minv && v < _maxv) {
      int i = int((v-_minv)/_dv); 
      if (i >= 0 && i < int(_nv)) {return unsigned(i);} 
    }
    // otherwise...
    return _nv;
  }


  /// return the mean value of all events given to histogram
  /// including those that were outside the histogram edges
  //double mean() const {return _weight_v / total_weight();}

  /// return the total weight in the histogram (inefficient)...
  double total_weight() const {
    if (!_have_total) {
      _total_weight = 0.0;
      for (unsigned i = 0; i < _weights.size(); i++) {
        _total_weight += _weights[i];}
      _have_total = true;
    }
    return _total_weight;
  }

  void add_entry(double u, double v, double weight = 1.0) {
    //if (v >= _minv && v < _maxv) {
    //  int i = int((v-_minv)/_dv); 
    //  if (i >= 0 && i < int(_weights.size())) _weights[i] += weight;
    //}
    _have_total = false;
    _weights[getbin(u_bin(u), v_bin(v))] += weight;
    //_weight_v += weight * v;
    //_weight_vsq += weight * v * v;
  };

  // Operations with constants ---------------------------------------
  SimpleHist2D & operator*=(double fact) {
    for (unsigned i = 0; i < outflow_size(); i++) (*this)[i] *= fact;
    //_weight_v *= fact;
    //_weight_vsq *= fact;
    _total_weight *= fact;
    return *this;
  };
  SimpleHist2D & operator/=(double fact) {
    *this *= 1.0/fact;
    return *this;
  };

  // Operations with another histogram -------------------------------
  SimpleHist2D & operator*=(const SimpleHist2D & other) {
    assert(other.outflow_size() == outflow_size());
    for (unsigned i = 0; i < outflow_size(); i++) (*this)[i] *= other[i];
    return *this;
  };

  SimpleHist2D & operator/=(const SimpleHist2D & other) {
    assert(other.outflow_size() == outflow_size());
    for (unsigned i = 0; i < outflow_size(); i++) (*this)[i] /= other[i];
    return *this;
  };

  SimpleHist2D & operator+=(const SimpleHist2D & other) {
    assert(other.outflow_size() == outflow_size());
    for (unsigned i = 0; i < outflow_size(); i++) (*this)[i] += other[i];
    //_weight_v += other._weight_v;
    //_weight_vsq += other._weight_vsq;
    if (_have_total && other._have_total) {
      _total_weight += other._total_weight;
    } else {_have_total = false;}
    return *this;
  };

  SimpleHist2D & operator-=(const SimpleHist2D & other) {
    assert(other.outflow_size() == outflow_size());
    for (unsigned i = 0; i < outflow_size(); i++) (*this)[i] -= other[i];
    return *this;
  };


  friend SimpleHist2D operator*(const SimpleHist2D & hist, double fact);
  friend SimpleHist2D operator/(const SimpleHist2D & hist, double fact);

private:
  double _minv, _maxv, _dv;
  double _minu, _maxu, _du;
  unsigned int _nu, _nv;
  std::valarray<double> _weights;
  std::string _name;
  //double _weight_v, _weight_vsq;
  mutable double _total_weight;
  mutable bool   _have_total;
};



// Binary operations with constants -----------------------------
inline SimpleHist2D operator*(const SimpleHist2D & hist, double fact) {
  SimpleHist2D result(hist.u_min(), hist.u_max(), hist.nu(),
		      hist.v_min(), hist.v_max(), hist.nv());
  for (unsigned i = 0; i < hist.outflow_size(); i++) result[i] = hist[i] * fact;
  //result._weight_v = hist._weight_v * fact;
  //result._weight_vsq = hist._weight_vsq * fact;
  return result;
}

inline SimpleHist2D operator/(const SimpleHist2D & hist, double fact) {
  SimpleHist2D result(hist.u_min(), hist.u_max(), hist.nu(),
		      hist.v_min(), hist.v_max(), hist.nv());
  for (unsigned i = 0; i < hist.outflow_size(); i++) result[i] = hist[i] / fact;
  //result._weight_v = hist._weight_v / fact;
  //result._weight_vsq = hist._weight_vsq / fact;
  return result;
}

inline SimpleHist2D operator*(double fact, const SimpleHist2D & hist) {
  return hist*fact;
}

inline SimpleHist2D operator/(double fact, const SimpleHist2D & hist) {
  return hist/fact;
}


// Binary operations with other histograms ------------------------
inline SimpleHist2D operator*(const SimpleHist2D & hista, const SimpleHist2D & histb) {
  assert(hista.outflow_size() == histb.outflow_size());
  SimpleHist2D result(hista.u_min(), hista.u_max(), hista.nu(),
		      hista.v_min(), hista.v_max(), hista.nv());
  for (unsigned i = 0; i < hista.outflow_size(); i++) result[i] = hista[i] * histb[i];
  return result;
}
inline SimpleHist2D operator/(const SimpleHist2D & hista, const SimpleHist2D & histb) {
  assert(hista.outflow_size() == histb.outflow_size());
  SimpleHist2D result(hista.u_min(), hista.u_max(), hista.nu(),
		      hista.v_min(), hista.v_max(), hista.nv());
  for (unsigned i = 0; i < hista.outflow_size(); i++) result[i] = hista[i] / histb[i];
  return result;
}
inline SimpleHist2D operator+(const SimpleHist2D & hista, const SimpleHist2D & histb) {
  assert(hista.outflow_size() == histb.outflow_size());
  SimpleHist2D result(hista.u_min(), hista.u_max(), hista.nu(),
		      hista.v_min(), hista.v_max(), hista.nv());
  for (unsigned i = 0; i < hista.outflow_size(); i++) result[i] = hista[i] + histb[i];
  return result;
}
inline SimpleHist2D operator-(const SimpleHist2D & hista, const SimpleHist2D & histb) {
  assert(hista.outflow_size() == histb.outflow_size());
  SimpleHist2D result(hista.u_min(), hista.u_max(), hista.nu(),
		      hista.v_min(), hista.v_max(), hista.nv());
  for (unsigned i = 0; i < hista.outflow_size(); i++) result[i] = hista[i] - histb[i];
  return result;
}


// Unary mathematical functions
inline SimpleHist2D sqrt(const SimpleHist2D & hist) {
  SimpleHist2D result(hist.u_min(), hist.u_max(), hist.nu(),
		      hist.v_min(), hist.v_max(), hist.nv());
  for (unsigned i = 0; i < hist.outflow_size(); i++) result[i] = sqrt(hist[i]);
  return result;
}

// Unary mathematical functions
inline SimpleHist2D pow2(const SimpleHist2D & hist) {
  SimpleHist2D result(hist.u_min(), hist.u_max(), hist.nu(),
		      hist.v_min(), hist.v_max(), hist.nv());
  for (unsigned i = 0; i < hist.outflow_size(); i++) result[i] = hist[i]*hist[i];
  return result;
}

/// output the histogram to standard output -- an operator<< might
/// have seemed nice, but less easy to generalize to multiple
/// histograms; the output is multipled by the factor norm.
inline void output(const SimpleHist2D & hist0, 
                   std::ostream * ostr = (&std::cout),
                   double norm = 1.0) {
  for (unsigned iu = 0; iu < hist0.nu(); iu++) {
  for (unsigned iv = 0; iv < hist0.nv(); iv++) {
    *ostr << hist0.u_binlo(iu)  << " " 
          << hist0.u_binmid(iu) << " "
          << hist0.u_binhi(iu) << " "
	  << hist0.v_binlo(iv)  << " " 
          << hist0.v_binmid(iv) << " "
          << hist0.v_binhi(iv) << " "
          << hist0(iu,iv)*norm << std::endl;
  }
  // make it readable by gnuplot
  *ostr << std::endl;
  }
}

inline void output(const SimpleHist2D & hist0, 
                   const SimpleHist2D & hist1, 
                   std::ostream * ostr = (&std::cout),
                   double norm = 1.0) {
  // assert(hist0.nv() == hist1.nv() && 
  //        hist0.u_min()  == hist1.u_min() &&
  //        hist0.u_max()  == hist1.u_max());
  for (unsigned iu = 0; iu < hist0.nu(); iu++) {
  for (unsigned iv = 0; iv < hist0.nv(); iv++) {
    *ostr << hist0.u_binlo(iu)  << " " 
          << hist0.u_binmid(iu) << " "
          << hist0.u_binhi(iu) << " "
	  << hist0.v_binlo(iv)  << " " 
          << hist0.v_binmid(iv) << " "
          << hist0.v_binhi(iv) << " "
          << hist0(iu,iv)*norm  << " "
          << hist1(iu,iv)*norm 
          << std::endl;
  }
  // make it readable by gnuplot
  *ostr << std::endl;
  }
}

inline void output(const SimpleHist2D & hist0, 
                   const SimpleHist2D & hist1, 
                   const SimpleHist2D & hist2, 
                   std::ostream * ostr = (&std::cout),
                   double norm = 1.0) {
  // assert(hist0.nv() == hist1.nv() && 
  //        hist0.u_min()  == hist1.u_min() &&
  //        hist0.u_max()  == hist1.u_max());
  for (unsigned iu = 0; iu < hist0.nu(); iu++) {
  for (unsigned iv = 0; iv < hist0.nv(); iv++) {
    *ostr << hist0.u_binlo(iu)  << " " 
          << hist0.u_binmid(iu) << " "
          << hist0.u_binhi(iu) << " "
	  << hist0.v_binlo(iv)  << " " 
          << hist0.v_binmid(iv) << " "
          << hist0.v_binhi(iv) << " "
          << hist0(iu,iv)*norm  << " "
          << hist1(iu,iv)*norm  << " "
          << hist2(iu,iv)*norm 
          << std::endl;
  }
  // make it readable by gnuplot
  *ostr << std::endl;
  }
}

inline void output(const SimpleHist2D & hist0, 
                   const SimpleHist2D & hist1, 
                   const SimpleHist2D & hist2, 
                   const SimpleHist2D & hist3, 
                   std::ostream * ostr = (&std::cout),
                   double norm = 1.0) {
  // assert(hist0.nv() == hist1.nv() && 
  //        hist0.u_min()  == hist1.u_min() &&
  //        hist0.u_max()  == hist1.u_max());
  for (unsigned iu = 0; iu < hist0.nu(); iu++) {
  for (unsigned iv = 0; iv < hist0.nv(); iv++) {
    *ostr << hist0.u_binlo(iu)  << " " 
          << hist0.u_binmid(iu) << " "
          << hist0.u_binhi(iu) << " "
	  << hist0.v_binlo(iv)  << " " 
          << hist0.v_binmid(iv) << " "
          << hist0.v_binhi(iv) << " "
          << hist0(iu,iv)*norm  << " "
          << hist1(iu,iv)*norm  << " "
          << hist2(iu,iv)*norm  << " "
          << hist3(iu,iv)*norm 
          << std::endl;
  }
  // make it readable by gnuplot
  *ostr << std::endl;
  }
}


// inline void output(const SimpleHist2D & hist0, 
//                    const SimpleHist2D & hist1, 
//                    std::ostream * ostr = (&std::cout),
//                    double norm = 1.0) {
//   assert(hist0.size() == hist1.size() && 
//          hist0.min()  == hist1.min() &&
//          hist0.max()  == hist1.max());
//   for (unsigned i = 0; i < hist0.size(); i++) {
//     *ostr << hist0.binlo(i)  << " " 
//           << hist0.binmid(i) << " "
//           << hist0.binhi(i) << " "
//           << hist0[i] * norm << " "
//           << hist1[i] * norm << " "
//           << std::endl;
//   }
// }
// 
// inline void output(const SimpleHist2D & hist0, 
//                    const SimpleHist2D & hist1, 
//                    const SimpleHist2D & hist2, 
//                    std::ostream * ostr = (&std::cout),
//                    double norm = 1.0) {
//   assert(hist0.size() == hist1.size() && 
//          hist0.min()  == hist1.min() &&
//          hist0.max()  == hist1.max());
//   assert(hist0.size() == hist2.size() && 
//          hist0.min()  == hist2.min() &&
//          hist0.max()  == hist2.max());
//   for (unsigned i = 0; i < hist0.size(); i++) {
//     *ostr << hist0.binlo(i)  << " " 
//           << hist0.binmid(i) << " "
//           << hist0.binhi(i) << " "
//           << hist0[i] * norm << " "
//           << hist1[i] * norm << " "
//           << hist2[i] * norm << " "
//           << std::endl;
//   }
// }
// 
// 
// inline void output(const SimpleHist2D & hist0, 
//                    const SimpleHist2D & hist1, 
//                    const SimpleHist2D & hist2, 
//                    const SimpleHist2D & hist3, 
//                    std::ostream * ostr = (&std::cout),
//                    double norm = 1.0) {
//   assert(hist0.size() == hist1.size() && 
//          hist0.min()  == hist1.min() &&
//          hist0.max()  == hist1.max());
//   assert(hist0.size() == hist2.size() && 
//          hist0.min()  == hist2.min() &&
//          hist0.max()  == hist2.max());
//   assert(hist0.size() == hist3.size() && 
//          hist0.min()  == hist3.min() &&
//          hist0.max()  == hist3.max());
//   for (unsigned i = 0; i < hist0.size(); i++) {
//     *ostr << hist0.binlo(i)  << " " 
//           << hist0.binmid(i) << " "
//           << hist0.binhi(i) << " "
//           << hist0[i] * norm << " "
//           << hist1[i] * norm << " "
//           << hist2[i] * norm << " "
//           << hist3[i] * norm << " "
//           << std::endl;
//   }
// }
// 
// inline void output(const SimpleHist2D & hist0, 
//                    const SimpleHist2D & hist1, 
//                    const SimpleHist2D & hist2, 
//                    const SimpleHist2D & hist3, 
//                    const SimpleHist2D & hist4, 
//                    std::ostream * ostr = (&std::cout),
//                    double norm = 1.0) {
//   assert(hist0.size() == hist1.size() && 
//          hist0.min()  == hist1.min() &&
//          hist0.max()  == hist1.max());
//   assert(hist0.size() == hist2.size() && 
//          hist0.min()  == hist2.min() &&
//          hist0.max()  == hist2.max());
//   assert(hist0.size() == hist3.size() && 
//          hist0.min()  == hist3.min() &&
//          hist0.max()  == hist3.max());
//   assert(hist0.size() == hist4.size() && 
//          hist0.min()  == hist4.min() &&
//          hist0.max()  == hist4.max());
//   for (unsigned i = 0; i < hist0.size(); i++) {
//     *ostr << hist0.binlo(i)  << " " 
//           << hist0.binmid(i) << " "
//           << hist0.binhi(i) << " "
//           << hist0[i] * norm << " "
//           << hist1[i] * norm << " "
//           << hist2[i] * norm << " "
//           << hist3[i] * norm << " "
//           << hist4[i] * norm << " "
//           << std::endl;
//   }
// }
// 
// inline void output(const SimpleHist2D & hist0, 
//                    const SimpleHist2D & hist1, 
//                    const SimpleHist2D & hist2, 
//                    const SimpleHist2D & hist3, 
//                    const SimpleHist2D & hist4, 
//                    const SimpleHist2D & hist5, 
//                    std::ostream * ostr = (&std::cout),
//                    double norm = 1.0) {
//   assert(hist0.size() == hist1.size() && 
//          hist0.min()  == hist1.min() &&
//          hist0.max()  == hist1.max());
//   assert(hist0.size() == hist2.size() && 
//          hist0.min()  == hist2.min() &&
//          hist0.max()  == hist2.max());
//   assert(hist0.size() == hist3.size() && 
//          hist0.min()  == hist3.min() &&
//          hist0.max()  == hist3.max());
//   assert(hist0.size() == hist4.size() && 
//          hist0.min()  == hist4.min() &&
//          hist0.max()  == hist4.max());
//   assert(hist0.size() == hist5.size() && 
//          hist0.min()  == hist5.min() &&
//          hist0.max()  == hist5.max());
//   for (unsigned i = 0; i < hist0.size(); i++) {
//     *ostr << hist0.binlo(i)  << " " 
//           << hist0.binmid(i) << " "
//           << hist0.binhi(i)  << " "
//           << hist0[i] * norm << " "
//           << hist1[i] * norm << " "
//           << hist2[i] * norm << " "
//           << hist3[i] * norm << " "
//           << hist4[i] * norm << " "
//           << hist5[i] * norm << " "
//           << std::endl;
//   }
// }
// 
// inline void output(const SimpleHist2D *hists, 
//                    const unsigned int nb_hist, 
//                    std::ostream * ostr = (&std::cout),
//                    double norm = 1.0) {
//   unsigned int ih;
//   for (ih=1;ih<nb_hist;ih++)
//     assert(hists[0].size() == hists[ih].size() && 
// 	   hists[0].min()  == hists[ih].min() &&
// 	   hists[0].max()  == hists[ih].max());
// 
//   for (unsigned i = 0; i < hists[0].size(); i++) {
//     *ostr << hists[0].binlo(i)  << " " 
//           << hists[0].binmid(i) << " "
//           << hists[0].binhi(i)  ;
//     for (ih=0;ih<nb_hist;ih++)
//       *ostr << " " << hists[ih][i] * norm;
//     *ostr << std::endl;
//   }
// }

#endif // __SIMPLEHIST2D_HH__
