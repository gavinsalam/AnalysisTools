# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]
### Added
    - a set of unit-tests for `SimpleHist` in unit-test/
    - new output routines for `SimpleHist`: `output(...)`
      `output_cumul(...)`, `output_diff(...)`, and
      `output_header(...)`
      
### Changed
    - SimpleHist now has underflow and overflow bins (access indices
      with `underflow_bin()` and `overflow_bin()` and the bins
      themselves as `underflow()` and `overflow()`
    - `SimpleHist::outflow()` no longer has the non-const reference option


## [0.1.0] - 2019-03-13
   An arbitrary starting point for this ChangeLog
   

