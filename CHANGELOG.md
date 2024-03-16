# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.2.0 devel]
- added view-hfile.py
- hfile.py has new features
  * HFile class for containing a whole file
  * Histogram class that contains a histogram with some attempt
    at automatically handling the various kinds of histograms
    produced by the Analysis Tools and for plotting
    (matplotlib should be loaded only if needed)
  * line_and_band(...) for plotting a line and a band

- for TemplateDefaultHist, added a declare_once option that returns the underlying
  histogram, e.g. `hists_err["name"].declare_once(binning).add_entry(value);`
- AveragingHist output now includes a "cols: ..." line
- added -git-info yes/no (or -no-git-info) option to determine whether
  runtime git info is included in the output
  
## [1.1.0] updates on 2023-01-31

### Changed

- automatic output interval now limited in lower and upper time
  between outputs through AnalysisBase::_min_time_between_writes and _max_time_between_writes

- combine-runs.pl now handles numbers that follow " +- " and " ± "  as errors
  (this is labeled as release 2.3.0 of combine-runs.pl)

- new virtual AnalysisBase::units_string() replaces earlier _units_string and 

- updated Catch to 3.1.0

- added -dump-argfile option

## [Unreleased]  2023-03-14
### Added
- many utility scripts (such as combine-runs.pl, hfile.py, etc.)
  are now included in the repository

## [Unreleased]
### Added
    - added hists_2d_compact, which output only bin centres
    - a set of unit-tests for `SimpleHist` in unit-test/
    - new output routines for `SimpleHist`: `output(...)`
      `output_cumul(...)`, `output_diff(...)`, and
      `output_header(...)`
    - GSLRandom has options for i/o of generator state
    
    
### Changed
    - SimpleHist now has underflow and overflow bins (access indices
      with `underflow_bin()` and `overflow_bin()` and the bins
      themselves as `underflow()` and `overflow()`
    - `SimpleHist::outflow()` no longer has the non-const reference option


## [0.1.0] - 2019-03-13
   An arbitrary starting point for this ChangeLog
   

