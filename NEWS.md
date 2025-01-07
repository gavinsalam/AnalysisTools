
# Towards 1.2.0

## New features in C++
- histogram objects in AnalysisBase now have a declare_once() method,
  which allows for usage such as `hists["some_name"].declare_once(binning).add_entry(value, weight);`
- AnalysisBase provides a `-git-info [yes|no]` option to control whether git info is
  determined at run time (default on, turn it off on systems where git would be slow)

## New features in Python (still in development)
- hfile.Histogram objects, which accommodate various kinds of histogram and know how to plot themselves
- hfile.HFile to read in a whole file, providing a list of hfile.Histogram objects
- new python/view-hfile.py script, relying on the above, to plot all histograms from a file


## Bug fixes
- averages[...].set_ref(xsc_label) was being ignored if the average was being
  filled with add_entry(value,weight); that's now fixed

# 1.1.0 2024-01-31

NEWS is incomplete up to and including 1.1.0
