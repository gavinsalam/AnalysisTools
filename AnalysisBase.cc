#include "AnalysisBase.hh"

template<> double DefaultHist::_lo =   0.0;
template<> double DefaultHist::_hi = 100.0;
template<> double DefaultHist::_bin_size = 1.0;

template<> double DefaultAveragingHist::_lo =   0.0;
template<> double DefaultAveragingHist::_hi = 100.0;
template<> double DefaultAveragingHist::_bin_size = 1.0;

template<> double DefaultCorrelationHist::_lo =   0.0;
template<> double DefaultCorrelationHist::_hi = 100.0;
template<> double DefaultCorrelationHist::_bin_size = 1.0;

