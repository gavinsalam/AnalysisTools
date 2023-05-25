#ifndef __SIMPLEHIST_HH__
#define __SIMPLEHIST_HH__

//STARTHEADER
// $Id: SimpleHist.hh 3392 2012-08-01 10:42:51Z salam $
//
// Copyright (c) 2007-2011, Matteo Cacciari, Gavin Salam and Gregory Soyez
//
//----------------------------------------------------------------------
//
//  SimpleHist is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  SimpleHist is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with FastJet; if not, write to the Free Software
//  Foundation, Inc.:
//      59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
//----------------------------------------------------------------------
//ENDHEADER


#include<valarray>
#include<string>
#include<cmath>
#include<iostream>
#include<vector>
#include "Binning.hh"

class SimpleHist : public Binning {
public:
  SimpleHist() {}
  SimpleHist(const Binning & binning): Binning(binning) {_init();}
  SimpleHist(double minv, double maxv, unsigned int n) : Binning (minv,maxv,n) {_init();}

  // NB: the int case is needed because if one passes an
  // int then compiler doesn't know whether to convert
  // to unsigned or double.
  SimpleHist(double minv, double maxv, int n) : Binning (minv,maxv,n) {
    _init();
  }

  SimpleHist(double minv, double maxv, double bin_size) : Binning (minv,maxv,bin_size)  {
    _init();
  }

  /// initialise the histogram, assuming the binning has already been set up
  void _init() override {
    _weights.resize(outflow_size());
    reset();
  }
  
  /// reset the contents of the histogram to zero (does not
  /// modify the histogram bounds)
  void reset() {
    _weights = 0.0;
    _weight_v = 0.0;
    _weight_vsq = 0.0;
    _total_weight = 0.0;
    _have_total = false;
    _n_entries = 0.0;
  }

  double & operator[](int i) {_have_total = false; return _weights[i];};
  const double & operator[](int i) const {return _weights[i];};

  /// returns the outflow bins
  double & underflow() {return _weights[underflow_bin()];}
  double & overflow()  {return _weights[overflow_bin()];}
  const double & underflow() const {return _weights[underflow_bin()];}
  const double & overflow()  const {return _weights[overflow_bin()];}
  
  /// returns the total outflow (underflow and overflow)
  double outflow() const {return underflow() + overflow();}

  /// return the mean value of all events given to histogram
  /// including those that were outside the histogram edges
  double mean() const {return _weight_v / total_weight();}

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

  double n_entries() const {return _n_entries;}

  /// add an entry to the bin in which v falls
  void add_entry(double v, double weight = 1.0) {
    _add_entry_ibin(v, bin(v),weight);
  };

  // Operations with constants ---------------------------------------
  SimpleHist & operator*=(double fact) {
    for (unsigned i = 0; i < outflow_size(); i++) (*this)[i] *= fact;
    _weight_v *= fact;
    _weight_vsq *= fact;
    _total_weight *= fact;
    return *this;
  }
  SimpleHist & operator/=(double fact) {
    *this *= 1.0/fact;
    return *this;
  }

  // Operations with another histogram -------------------------------
  SimpleHist & operator*=(const SimpleHist & other) {
    assert(other.outflow_size() == outflow_size());
    for (unsigned i = 0; i < outflow_size(); i++) (*this)[i] *= other[i];
    return *this;
  };

  SimpleHist & operator/=(const SimpleHist & other) {
    assert(other.outflow_size() == outflow_size());
    for (unsigned i = 0; i < outflow_size(); i++) (*this)[i] /= other[i];
    return *this;
  };

  SimpleHist & operator+=(const SimpleHist & other) {
    assert(other.outflow_size() == outflow_size());
    for (unsigned i = 0; i < outflow_size(); i++) (*this)[i] += other[i];
    _weight_v += other._weight_v;
    _weight_vsq += other._weight_vsq;
    if (_have_total && other._have_total) {
      _total_weight += other._total_weight;
    } else {_have_total = false;}
    return *this;
  };

  SimpleHist & operator-=(const SimpleHist & other) {
    assert(other.outflow_size() == outflow_size());
    for (unsigned i = 0; i < outflow_size(); i++) (*this)[i] -= other[i];
    return *this;
  };

  // output operations ----------------------------------------------
  /// output a header with the total weight and the mean value of the histogram
  std::ostream & output_total_and_mean(std::ostream & ostr, const std::string & prefix = "# ") {
    ostr << prefix << "total_weight = " << total_weight() << std::endl;
    ostr << prefix << "n_entries = " << n_entries() << std::endl;
    ostr << prefix << "mean = " << mean() << std::endl;
    return ostr;
  }

  /// output the histogram, including the header.
  /// The prefix is used for lines not intended to be used in normal plots
  std::ostream & output(std::ostream & ostr, const std::string & prefix = "# ") {
    output_total_and_mean(ostr, prefix);
    ostr << prefix << "cols: vlo vmid vhi hist[v]" << std::endl;
    ostr << prefix << "under flow bin " << underflow() << std::endl;
    for (unsigned i = 0; i < size() ; i++) {
      ostr << binlo(i)  << " " 
           << binmid(i) << " "
           << binhi(i) << " "
           << (*this)[i] << std::endl;
    }
    ostr << prefix << "over flow bin " << overflow() << std::endl;
    return ostr;
  }
  
  /// output the histogram, including the header, as a differential
  /// histogram d/dv. Outflow bins don't get any bin width normalisaion
  /// The prefix is used for lines not intended to be used in normal plots
  std::ostream & output_diff(std::ostream & ostr, const std::string & prefix = "# ") {
    output_total_and_mean(ostr, prefix);
    ostr << prefix << "cols: vlo vmid vhi dhist/dv (outflow not normalised)" << std::endl;
    ostr << prefix << "under flow bin " << underflow() << std::endl;
    for (unsigned i = 0; i < size() ; i++) {
      ostr << binlo(i)  << " " 
           << binmid(i) << " "
           << binhi(i) << " "
           << (*this)[i] / (binhi(i) - binlo(i)) << std::endl;
    }
    ostr << prefix << "over flow bin " << overflow() << std::endl;
    return ostr;
  }
  
  /// output the cumulative histogram.
  /// The prefix is used for lines not intended to be used in normal plots
  std::ostream & output_cumul(std::ostream & ostr, const std::string & prefix = "# ") {
    output_total_and_mean(ostr, prefix);
    double cumul = underflow();
    ostr << prefix << "cols: v hist_integral_up_to_v" << std::endl;
    ostr << min() << " " << cumul << std::endl;
    for (unsigned i = 0; i < size(); i++) {
      cumul += (*this)[i];
      ostr << binhi(i) << " " << cumul << std::endl;
    }
    cumul += overflow();
    ostr << prefix << "with_overflow " << cumul << std::endl;
    return ostr;
  }
  
  friend SimpleHist operator*(const SimpleHist & hist, double fact);
  friend SimpleHist operator/(const SimpleHist & hist, double fact);

protected:

  /// add an entry with the bin already worked out
  /// (NB: this is virtual so that the _add_entry_ibin call in
  /// add_entry can automatically call derived class add_entry_ibin
  /// functions).  IT DOES NOT CHECK THAT v and i are consistent
  virtual void _add_entry_ibin(double v, unsigned ibin, double weight) {
    _have_total = false;
    _weights[ibin] += weight;
    _weight_v += weight * v;
    _weight_vsq += weight * v * v;
    _n_entries += 1.0;
  }

  
  std::valarray<double> _weights;
  std::string _name;
  double _weight_v, _weight_vsq;
  mutable double _total_weight;
  mutable bool   _have_total;
  double _n_entries;
};



// Binary operations with constants -----------------------------
inline SimpleHist operator*(const SimpleHist & hist, double fact) {
  SimpleHist result = hist;
  result *= fact;
  return result;

  //SimpleHist result(hist.min(), hist.max(), hist.size());
  //for (unsigned i = 0; i < hist.outflow_size(); i++) result[i] = hist[i] * fact;
  //result._weight_v = hist._weight_v * fact;
  //result._weight_vsq = hist._weight_vsq * fact;
  //result._have_total = hist.have_total;
  //return result;
}

inline SimpleHist operator/(const SimpleHist & hist, double fact) {
  SimpleHist result = hist;
  result *= fact;
  return result;

  //SimpleHist result(hist.min(), hist.max(), hist.size());
  //for (unsigned i = 0; i < hist.outflow_size(); i++) result[i] = hist[i] / fact;
  //result._weight_v = hist._weight_v / fact;
  //result._weight_vsq = hist._weight_vsq / fact;
  //return result;
}

inline SimpleHist operator*(double fact, const SimpleHist & hist) {
  return hist*fact;
}

inline SimpleHist operator/(double fact, const SimpleHist & hist) {
  return hist/fact;
}


// Binary operations with other histograms ------------------------
inline SimpleHist operator*(const SimpleHist & hista, const SimpleHist & histb) {
  assert(hista.outflow_size() == histb.outflow_size());
  SimpleHist result(hista.min(), hista.max(), hista.size());
  for (unsigned i = 0; i < hista.outflow_size(); i++) result[i] = hista[i] * histb[i];
  return result;
}
inline SimpleHist operator/(const SimpleHist & hista, const SimpleHist & histb) {
  assert(hista.outflow_size() == histb.outflow_size());
  SimpleHist result(hista.min(), hista.max(), hista.size());
  for (unsigned i = 0; i < hista.outflow_size(); i++) result[i] = hista[i] / histb[i];
  return result;
}
inline SimpleHist operator+(const SimpleHist & hista, const SimpleHist & histb) {
  assert(hista.outflow_size() == histb.outflow_size());
  SimpleHist result(hista.min(), hista.max(), hista.size());
  for (unsigned i = 0; i < hista.outflow_size(); i++) result[i] = hista[i] + histb[i];
  return result;
}
inline SimpleHist operator-(const SimpleHist & hista, const SimpleHist & histb) {
  assert(hista.outflow_size() == histb.outflow_size());
  SimpleHist result(hista.min(), hista.max(), hista.size());
  for (unsigned i = 0; i < hista.outflow_size(); i++) result[i] = hista[i] - histb[i];
  return result;
}


// Unary mathematical functions
inline SimpleHist sqrt(const SimpleHist & hist) {
  SimpleHist result(hist.min(), hist.max(), hist.size());
  for (unsigned i = 0; i < hist.outflow_size(); i++) result[i] = sqrt(hist[i]);
  return result;
}

// Unary mathematical functions
inline SimpleHist pow2(const SimpleHist & hist) {
  SimpleHist result(hist.min(), hist.max(), hist.size());
  for (unsigned i = 0; i < hist.outflow_size(); i++) result[i] = hist[i]*hist[i];
  return result;
}

/// output the histogram to standard output -- an operator<< might
/// have seemed nice, but less easy to generalize to multiple
/// histograms; the output is multipled by the factor norm.
inline void output(const SimpleHist & hist0, 
                   std::ostream * ostr = (&std::cout),
                   double norm = 1.0) {
  for (unsigned i = 0; i < hist0.size(); i++) {
    *ostr << hist0.binlo(i)  << " " 
          << hist0.binmid(i) << " "
          << hist0.binhi(i) << " "
          << hist0[i]*norm << std::endl;
  }
}


inline std::ostream & operator<<(std::ostream & ostr, const SimpleHist & hist) {
  output(hist, &ostr);
  return ostr;
}

inline void output(const SimpleHist & hist0, 
                   const SimpleHist & hist1, 
                   std::ostream * ostr = (&std::cout),
                   double norm = 1.0) {
  assert(hist0.size() == hist1.size() && 
         hist0.min()  == hist1.min() &&
         hist0.max()  == hist1.max());
  for (unsigned i = 0; i < hist0.size(); i++) {
    *ostr << hist0.binlo(i)  << " " 
          << hist0.binmid(i) << " "
          << hist0.binhi(i) << " "
          << hist0[i] * norm << " "
          << hist1[i] * norm << " "
          << std::endl;
  }
}

inline void output(const SimpleHist & hist0, 
                   const SimpleHist & hist1, 
                   const SimpleHist & hist2, 
                   std::ostream * ostr = (&std::cout),
                   double norm = 1.0) {
  assert(hist0.size() == hist1.size() && 
         hist0.min()  == hist1.min() &&
         hist0.max()  == hist1.max());
  assert(hist0.size() == hist2.size() && 
         hist0.min()  == hist2.min() &&
         hist0.max()  == hist2.max());
  for (unsigned i = 0; i < hist0.size(); i++) {
    *ostr << hist0.binlo(i)  << " " 
          << hist0.binmid(i) << " "
          << hist0.binhi(i) << " "
          << hist0[i] * norm << " "
          << hist1[i] * norm << " "
          << hist2[i] * norm << " "
          << std::endl;
  }
}


inline void output(const SimpleHist & hist0, 
                   const SimpleHist & hist1, 
                   const SimpleHist & hist2, 
                   const SimpleHist & hist3, 
                   std::ostream * ostr = (&std::cout),
                   double norm = 1.0) {
  assert(hist0.size() == hist1.size() && 
         hist0.min()  == hist1.min() &&
         hist0.max()  == hist1.max());
  assert(hist0.size() == hist2.size() && 
         hist0.min()  == hist2.min() &&
         hist0.max()  == hist2.max());
  assert(hist0.size() == hist3.size() && 
         hist0.min()  == hist3.min() &&
         hist0.max()  == hist3.max());
  for (unsigned i = 0; i < hist0.size(); i++) {
    *ostr << hist0.binlo(i)  << " " 
          << hist0.binmid(i) << " "
          << hist0.binhi(i) << " "
          << hist0[i] * norm << " "
          << hist1[i] * norm << " "
          << hist2[i] * norm << " "
          << hist3[i] * norm << " "
          << std::endl;
  }
}

inline void output(const SimpleHist & hist0, 
                   const SimpleHist & hist1, 
                   const SimpleHist & hist2, 
                   const SimpleHist & hist3, 
                   const SimpleHist & hist4, 
                   std::ostream * ostr = (&std::cout),
                   double norm = 1.0) {
  assert(hist0.size() == hist1.size() && 
         hist0.min()  == hist1.min() &&
         hist0.max()  == hist1.max());
  assert(hist0.size() == hist2.size() && 
         hist0.min()  == hist2.min() &&
         hist0.max()  == hist2.max());
  assert(hist0.size() == hist3.size() && 
         hist0.min()  == hist3.min() &&
         hist0.max()  == hist3.max());
  assert(hist0.size() == hist4.size() && 
         hist0.min()  == hist4.min() &&
         hist0.max()  == hist4.max());
  for (unsigned i = 0; i < hist0.size(); i++) {
    *ostr << hist0.binlo(i)  << " " 
          << hist0.binmid(i) << " "
          << hist0.binhi(i) << " "
          << hist0[i] * norm << " "
          << hist1[i] * norm << " "
          << hist2[i] * norm << " "
          << hist3[i] * norm << " "
          << hist4[i] * norm << " "
          << std::endl;
  }
}

inline void output(const SimpleHist & hist0, 
                   const SimpleHist & hist1, 
                   const SimpleHist & hist2, 
                   const SimpleHist & hist3, 
                   const SimpleHist & hist4, 
                   const SimpleHist & hist5, 
                   std::ostream * ostr = (&std::cout),
                   double norm = 1.0) {
  assert(hist0.size() == hist1.size() && 
         hist0.min()  == hist1.min() &&
         hist0.max()  == hist1.max());
  assert(hist0.size() == hist2.size() && 
         hist0.min()  == hist2.min() &&
         hist0.max()  == hist2.max());
  assert(hist0.size() == hist3.size() && 
         hist0.min()  == hist3.min() &&
         hist0.max()  == hist3.max());
  assert(hist0.size() == hist4.size() && 
         hist0.min()  == hist4.min() &&
         hist0.max()  == hist4.max());
  assert(hist0.size() == hist5.size() && 
         hist0.min()  == hist5.min() &&
         hist0.max()  == hist5.max());
  for (unsigned i = 0; i < hist0.size(); i++) {
    *ostr << hist0.binlo(i)  << " " 
          << hist0.binmid(i) << " "
          << hist0.binhi(i)  << " "
          << hist0[i] * norm << " "
          << hist1[i] * norm << " "
          << hist2[i] * norm << " "
          << hist3[i] * norm << " "
          << hist4[i] * norm << " "
          << hist5[i] * norm << " "
          << std::endl;
  }
}

inline void output(const SimpleHist *hists, 
                   const unsigned int nb_hist, 
                   std::ostream * ostr = (&std::cout),
                   double norm = 1.0) {
  unsigned int ih;
  for (ih=1;ih<nb_hist;ih++)
    assert(hists[0].size() == hists[ih].size() && 
	   hists[0].min()  == hists[ih].min() &&
	   hists[0].max()  == hists[ih].max());

  for (unsigned i = 0; i < hists[0].size(); i++) {
    *ostr << hists[0].binlo(i)  << " " 
          << hists[0].binmid(i) << " "
          << hists[0].binhi(i)  ;
    for (ih=0;ih<nb_hist;ih++)
      *ostr << " " << hists[ih][i] * norm;
    *ostr << std::endl;
  }
}

inline void output(const std::vector<SimpleHist*> &hists, 
                   std::ostream * ostr = (&std::cout),
                   double norm = 1.0) {
  unsigned int ih;
  unsigned int nb_hist = hists.size();
  for (ih=1;ih<nb_hist;ih++)
    assert(hists[0]->size() == hists[ih]->size() && 
	   hists[0]->min()  == hists[ih]->min() &&
	   hists[0]->max()  == hists[ih]->max());

  for (unsigned i = 0; i < hists[0]->size(); i++) {
    *ostr << hists[0]->binlo(i)  << " " 
          << hists[0]->binmid(i) << " "
          << hists[0]->binhi(i)  ;
    for (ih=0;ih<nb_hist;ih++)
      *ostr << " " << (*hists[ih])[i] * norm;
    *ostr << std::endl;
  }
}

inline void output(const std::vector<SimpleHist> &hists, 
                   std::ostream * ostr = (&std::cout),
                   double norm = 1.0) {
  unsigned int ih;
  unsigned int nb_hist = hists.size();
  for (ih=1;ih<nb_hist;ih++)
    assert(hists[0].size() == hists[ih].size() && 
	   hists[0].min()  == hists[ih].min() &&
	   hists[0].max()  == hists[ih].max());

  for (unsigned i = 0; i < hists[0].size(); i++) {
    *ostr << hists[0].binlo(i)  << " " 
          << hists[0].binmid(i) << " "
          << hists[0].binhi(i)  ;
    for (ih=0;ih<nb_hist;ih++)
      *ostr << " " << hists[ih][i] * norm;
    *ostr << std::endl;
  }
}

#endif // __SIMPLEHIST_HH__
