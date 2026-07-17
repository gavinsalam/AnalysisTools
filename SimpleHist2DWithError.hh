#pragma once 
#include "SimpleHist2D.hh"

class SimpleHist2DWithError : public SimpleHist2D {
public:


  SimpleHist2DWithError() {};
  SimpleHist2DWithError(double minu, double maxu, int nu,
			double minv, double maxv, int nv) {
    declare(minu, maxu, unsigned(nu), minv, maxv, unsigned(nv));
  }

  SimpleHist2DWithError(double minu, double maxu, unsigned int nu,
	       double minv, double maxv, unsigned int nv) {
    declare(minu, maxu, nu, minv, maxv, nv);
  }

  SimpleHist2DWithError(double minu, double maxu, double bin_size_u,
	       double minv, double maxv, double bin_size_v) {
    SimpleHist2D::declare(minu, maxu, bin_size_u, minv, maxv, bin_size_v);
  }

  void declare(double minu, double maxu, unsigned int nu,
	             double minv, double maxv, unsigned int nv) override {
    SimpleHist2D::declare(minu, maxu, nu, minv, maxv, nv);
    _weights_sumsqr.resize(outflow_size());
    _weights_sumsqr = 0.0;
  } 

  /// returns the error on the bin's contents
  double error(int iu, int iv) const {
    return error(getbin(iu,iv));
  }
  double error(unsigned i) const {return _error_calc(_weights[i], _weights_sumsqr[i]);}


  double sumsqr(int iu, int iv) const {
    return _weights_sumsqr[getbin(iu,iv)];
  }
  double sumsqr(unsigned i) const {return _weights_sumsqr[i];}

  // Operations with constants ---------------------------------------
  SimpleHist2DWithError & operator*=(double fact) {
    double factsqr = fact*fact;
    for (unsigned i = 0; i < outflow_size(); i++) {
      (*this)[i] *= fact;
      _weights_sumsqr[i] *= factsqr;
    }
    _total_weight *= fact;
    return *this;
  };
  SimpleHist2DWithError & operator/=(double fact) {
    *this *= 1.0/fact;
    return *this;
  };
  
  
  SimpleHist2DWithError & operator+=(const SimpleHist2DWithError & other) {
    assert(other.outflow_size() == outflow_size());
    for (unsigned i = 0; i < outflow_size(); i++) {
      (*this)[i] += other[i];
      _weights_sumsqr[i] += other._weights_sumsqr[i];
    }
    _n_entries += other._n_entries;
    if (_have_total && other._have_total) {
      _total_weight += other._total_weight;
    } else {_have_total = false;}
    return *this;
  };

  SimpleHist2DWithError & operator-=(const SimpleHist2DWithError & other) {
    auto minus_other = other;
    minus_other *= -1.0;
    *this += minus_other;
    return *this;
  }
  
protected:
  void _add_entry_ibin(unsigned int ibin, double weight) override {
    SimpleHist2D::_add_entry_ibin(ibin, weight);
    _weights_sumsqr[ibin] += weight*weight;
  }  

  double _error_calc(double sum, double sumsq) const {
    return std::sqrt(std::abs(sumsq - sum*sum/n_entries()));
  }

  std::valarray<double> _weights_sumsqr;

};

inline SimpleHist2DWithError operator+(const SimpleHist2DWithError & hista, const SimpleHist2DWithError & histb) {
  auto result = hista;
  result += histb;
  return result;
}

inline SimpleHist2DWithError operator-(const SimpleHist2DWithError & hista, const SimpleHist2DWithError & histb) {
  auto result = hista;
  result -= histb;
  return result;
}
inline SimpleHist2DWithError operator*(const SimpleHist2DWithError & hist, double fact) {  
  SimpleHist2DWithError result(hist);
  result *= fact;
  return result;
}
inline SimpleHist2DWithError operator/(const SimpleHist2DWithError & hist, double fact) {  
  return hist*(1.0/fact);
}

inline SimpleHist2DWithError operator*(double fact, const SimpleHist2DWithError & hist) {
  return hist*fact;
}
inline SimpleHist2DWithError operator/(double fact, const SimpleHist2DWithError & hist) {
  return hist/fact;
}
