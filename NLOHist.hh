#ifndef __NLOHIST_HH__
#define __NLOHIST_HH__

#include "SimpleHist.hh"
#include<iostream>
#include<vector>
#include<string>
#include<fstream>
#include<iomanip>

/// Class for carrying out histogramming in NLO calculations
///
/// 
class NLOHist {
public:
  NLOHist() {
    declare("dummy",0.0, 1.0, 10);
  }

  /// initialise a histogram for NLO studies
  NLOHist(const std::string & name, double minv, double maxv, int n) {
    declare(name, minv, maxv, n);
  }

  /// a way to get the histogram produced outside of a constructor
  void declare(const std::string & name, double minv, double maxv, int n) {
    _name = name;
    // initialise _main, _sqr, _tmp
    _main.declare(minv, maxv, n);
    _sqr .declare(minv, maxv, n);
    _tmp .declare(minv, maxv, n);
    // our normalisation
    _nev = 0.0;
    _nev_bad = 0.0;
    _event_is_bad = false;
  }

  /// return the number of events
  double nev() const {return _nev;}
  double nev_bad() const {return _nev_bad;}

  /// add an entry to the histogram
  void add_entry(double v, double weight) {
    // test for NaN
    if (weight != weight) {
      _event_is_bad = true;
      return;
    }
    // otherwise, bin it
    int ibin = _tmp.bin(v);
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
  void write(std::ostream & ostr, int index) {
    // make sure we're not in the middle of an event
    assert(_in_use.size() == 0);
    ostr << "# " << _name << " (index = " << index << ")"<< std::endl;
    // calculate the errors
    for (unsigned ibin = 0; ibin < _tmp.outflow_size(); ibin++) {
      // this will be used for the error, and contains a factor
      // nev**2 inside the sqrt currently (will be removed when writing)
      _tmp[ibin] = sqrt(std::abs(_sqr[ibin] - _main[ibin]*_main[ibin]/_nev));
    }

    // finally output the histogram with errors
    ostr << std::setprecision(8); // get decent precision
    output(_main, _tmp, &ostr, 1.0/_nev/_tmp.binsize());

    // reset tmp
    _tmp *= 0.0;
  }

  /// For Event 2, we want to be able to pass a normalisation (the
  /// normalisation should NOT include the bin size)
  ///
  /// Here, weights for each event are just x-sections/nev, so if one
  /// is not normalising to a cross-section, one should make sure that
  /// the normalisation includes a factor _nev/nev_total. This is left
  /// to the user. In other words, if the histogram is normalised to
  /// give a 1/sigma_tot dsigma/dX, the normalisation should be
  /// 1/sigma_tot. If instead, we want dsigma/dX, the normalisation
  /// should be nev_total_expected/nev_so_far
  void write(std::ostream & ostr, double norm) {
    // make sure we're not in the middle of an event
    assert(_in_use.size() == 0);
    ostr << "# " << _name << std::endl;
    // calculate the errors
    for (unsigned ibin = 0; ibin < _tmp.outflow_size(); ibin++) {
      // this will be used for the error, and contains a factor
      // nev**2 inside the sqrt currently (will be removed when writing)
      _tmp[ibin] = sqrt(std::abs(_sqr[ibin] - _main[ibin]*_main[ibin]/_nev));
    }

    // finally output the histogram with errors
    ostr << std::setprecision(8); // get decent precision
    output(_main, _tmp, &ostr, norm/_tmp.binsize());

    // reset tmp
    _tmp *= 0.0;
  }

  /// return the name given to the histogram
  const std::string & name() const {return _name;}
  const double min() const {return _main.min();}
  const double max() const {return _main.max();}

private:
  
  double      _nev, _nev_bad;
  SimpleHist  _main, _sqr, _tmp;
  std::vector<int> _in_use;
  std::string _name;
  bool        _event_is_bad;
};
#endif // __NLOHIST_HH__

