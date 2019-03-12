#ifndef __NLOHISTGENERIC_HH__
#define __NLOHISTGENERIC_HH__

#include "SimpleHist.hh"
#include<vector>
#include<string>
#include<ostream>
#include<iomanip>
#include<cmath>

using namespace std;

/// Class for carrying out histogramming in NLO calculations
///
/// 
class NLOHistGeneric {
public:
  NLOHistGeneric() {
    std::vector<double> lims;
    lims.push_back(0.0);
    lims.push_back(1.0);
    declare("", lims);
    _explicit_init_done = false; // allow for subsequent re-initialisation
  }

  /// initialise a histogram for NLO studies
  NLOHistGeneric(std::vector<double> limits) {
    declare("", limits);
  }

  /// initialise a histogram for NLO studies
  NLOHistGeneric(const std::string & name, std::vector<double> limits) {
    declare(name, limits);
  }

  /// a way to get the histogram produced outside of a constructor
  void declare(std::vector<double> limits) {
    declare("", limits);
  }

  void declare(double vmin, double vmax, double dv) {
    declare("", vmin, vmax, dv);
  }
  
  void declare(const std::string & name, double vmin, double vmax, double dv) {
    int nbin = abs(int(0.5+ (vmax - vmin)/dv));
    double new_dv = (vmax - vmin)/nbin;
    vector<double> limits(nbin+1);
    for (int i = 0; i <= nbin; i++) {
      limits[i] = vmin + new_dv*i;
    }
    declare(name, limits);
  }
  
  /// a way to get the histogram produced outside of a constructor
  void declare(const std::string & name, std::vector<double> limits) {
    _name = name;

    _limits = limits;

    assert(_limits.size() >= 2U);
    
    // initialise _main, _sqr, _tmp
    _main.resize(_limits.size()-1); fill(_main.begin(), _main.end(), 0.0);
    _sqr .resize(_limits.size()-1); fill(_sqr.begin(), _sqr.end(), 0.0);
    _tmp .resize(_limits.size()-1); fill(_tmp.begin(), _tmp.end(), 0.0);

    // our normalisation
    _nev = 0.0;
    _nev_bad = 0.0;
    _event_is_bad = false;

    _explicit_init_done = true;
  }

  /// return the number of events
  double nev() const {return _nev;}
  double nev_bad() const {return _nev_bad;}

  /// check in which bin we are
  int get_bin(double v){
    if (v<_limits[0]) return -1;
    if (v>_limits[_limits.size()-1]) return -1; // as an error instead of limits.size();

    int ibin=0;
    while (v>_limits[ibin+1]) ibin++;

    return ibin;
  }

  // declare the histogram limits and add an entry
  void set_lims_add_entry(double vmin, double vmax, double dv, double v, double weight) {
    if (!_explicit_init_done) {
      declare(vmin, vmax, dv);
    } else {
      // front limit should be identical
      assert(_limits.front() == vmin);
      // back limit may have rounding errors...
      const double tolerance = 1e-7;
      assert(abs(_limits.back() - vmax)
             < tolerance * std::max( abs(vmax), abs(_limits.back()-_limits.front())) );
    }
    add_entry(v, weight);
  }
  
  // declare the histogram limits and add an entry
  void set_lims_add_entry(const vector<double> & bin_edges, double v, double weight) {
    if (!_explicit_init_done) {
      declare(bin_edges);
    } else {
      assert(_limits.size()  == bin_edges.size());
      assert(_limits.front() == bin_edges.front());
      assert(_limits.back()  == bin_edges.back());
    }
    add_entry(v, weight);
  }

  /// add an entry to the histogram
  void add_entry(double v, double weight) {
    // test for NaN
    if (weight != weight) {
      _event_is_bad = true;
      return;
    }
    // otherwise, bin it
    int ibin = get_bin(v);
    if (ibin<0) return;

    if (_tmp[ibin] == 0) _in_use.push_back(ibin);
    _tmp[ibin] += weight;
  }


  /// collect the information binned within the event
  /// and transfer it to the definitive store (_main, _sqr)
  void collate_event() {
    _nev += 1;
    for (unsigned i = 0; i < _in_use.size(); i++) {
      int ibin = _in_use[i];
      // only transfer info from whole event, if there
      // were no bad events in it at all
      if (!_event_is_bad) {
        _main[ibin] += _tmp[ibin];
        _sqr [ibin] += _tmp[ibin]*_tmp[ibin];
      }
      _tmp[ibin]   = 0.0;
    }
    _in_use.clear();

    // if things were bad, reset histograms
    if (_event_is_bad) {
      _event_is_bad = false;
      _nev_bad++;
    }

  }


  /// output the histogram with errors, in a format similar to
  /// Zoltan's
  ///
  /// NB: we here normalise to nev(), but suspect that Zoltan
  ///     instead normalises to nev()-nev_bad(); 
  ///
  ///     we tend to believe that our choice is more correct [though
  ///     the effect should be sub per-mille level, at worst, for
  ///     dipole phase-space]: the normalisation of the "good" bins,
  ///     should not be affected by the fact that in "difficult" bins
  ///     the finite-precision numerics gave a NaN rather than the
  ///     correct answer.
  ///
  void write(std::ostream & ostr, int hist_index) const {
    double norm = 1/_nev;
    ostr << "# " << _name << " (index = " << hist_index << ")" << std::endl;
    write_norm(ostr, norm);
    
    // // make sure we're not in the middle of an event
    // assert(_in_use.size() == 0);
    // ostr << "# " << _name << " (index = " << index << ")"<< std::endl;
    // // calculate the errors
    // for (unsigned ibin = 0; ibin < _limits.size()-1; ibin++) {
    //   // this will be used for the error, and contains a factor
    //   // nev**2 inside the sqrt currently (will be removed when writing)
    //   _tmp[ibin] = sqrt(std::abs(_sqr[ibin] - _main[ibin]*_main[ibin]/_nev));
    // }
    // 
    // // finally output the histogram with errors
    // ostr << std::setprecision(8); // get decent precision
    // for (unsigned ibin = 0; ibin < _limits.size()-1; ibin++) {
    //   double binsize = _limits[ibin+1] - _limits[ibin];
    //   ostr << _limits[ibin] << " " << 0.5*(_limits[ibin]+_limits[ibin+1]) << " " << _limits[ibin+1] << " "
    //        << _main[ibin]/_nev/binsize << " " << _tmp[ibin]/_nev/binsize << std::endl;
    //   _tmp[ibin] = 0.0;
    // }
  }

  /// write the output with a user-chosen multiplicative normalisation
  //// - bin width is still automatically divided out
  ///  - this version leaves out the label line
  void write_norm(std::ostream & ostr, double norm) const {
    // make sure we're not in the middle of an event
    assert(_in_use.size() == 0);
    // calculate the errors
    std::vector<double> tmp(_tmp.size(), 0.0);
    for (unsigned ibin = 0; ibin < _limits.size()-1; ibin++) {
      // this will be used for the error, and contains a factor
      // nev**2 inside the sqrt currently (will be removed when writing)
      tmp[ibin] = sqrt(std::abs(_sqr[ibin] - _main[ibin]*_main[ibin]/_nev));
    }

    // finally output the histogram with errors
    ostr << std::setprecision(8); // get decent precision
    for (unsigned ibin = 0; ibin < _limits.size()-1; ibin++) {
      double binsize = _limits[ibin+1] - _limits[ibin];
      ostr << _limits[ibin] << " " << 0.5*(_limits[ibin]+_limits[ibin+1]) << " " << _limits[ibin+1] << " "
	   << _main[ibin]*norm /binsize << " " << tmp[ibin]*norm /binsize << std::endl;
      tmp[ibin] = 0.0;
    }
  }

  
  /// return the name given to the histogram
  const std::string & name() const {return _name;}
  const double min() const {return _limits[0];}
  const double max() const {return _limits[_limits.size()-1];}

private:

  bool _explicit_init_done;
  
  double      _nev, _nev_bad;

  std::vector<double> _limits;
  std::vector<double> _main, _sqr, _tmp;

  std::vector<int> _in_use;
  std::string _name;
  bool        _event_is_bad;
};

#endif // __NLOHISTGENERIC_HH__

