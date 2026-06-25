AnalysisTools
=============

A set of C++ header files (and some accompanying python) to provide
various histogramming, ntuple, etc. facilities.

Usage
-----

See tests/test-all-types.cc for an example of usage 
(run ./mkmk to create a Makefile).

Using AnalysisTools from another CMake project
----------------------------------------------

The supported CMake integration path is to add this repository to your
project source tree, for example as a git submodule.

This repository expects:

- CMake 3.16 or later
- a C++20 compiler
- GSL to be installed and discoverable by CMake

A minimal parent `CMakeLists.txt` can look like this:

```cmake
cmake_minimum_required(VERSION 3.16)
project(MyProject LANGUAGES CXX)

add_subdirectory(external/AnalysisTools)

add_executable(my-analysis main.cc)
target_link_libraries(my-analysis PRIVATE AnalysisTools::AnalysisTools)
```

This gives your target:

- the AnalysisTools headers
- the `CmdLine` dependency
- the GSL dependency
- the required C++20 compile feature

If you do not want AnalysisTools to build its own `tests/` and
`unit-tests/` targets when included as a subdirectory, you can set
`ANALYSISTOOLS_BUILD_TESTS` and `ANALYSISTOOLS_BUILD_UNIT_TESTS` to
`OFF` before calling `add_subdirectory(...)`. That is optional.

In your C++ source, include the headers directly, for example

```cc
#include "AnalysisBase.hh"
#include "GSLRandom.hh"
```

At present, the CMake support is aimed at `add_subdirectory(...)`
usage inside another project. The repository does not yet provide a
documented standalone install-and-`find_package(AnalysisTools)` flow.

Authors
-------

This code was originally written by Matteo Cacciari, Gavin Salam and
Gregory Soyez. It has since received contributions also from Keith 
Hamilton and Pier Monni.
