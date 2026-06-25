cmake_minimum_required(VERSION 3.16)

project(AnalysisTools VERSION 1.2.1 LANGUAGES CXX)

include(CTest)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(ANALYSISTOOLS_BUILD_TESTS "Build example programs in tests/" ON)
option(ANALYSISTOOLS_BUILD_UNIT_TESTS "Build unit tests in unit-tests/" ON)

set(CMDLINE_BUILD_EXAMPLE OFF CACHE BOOL "Build the CmdLine example program" FORCE)
set(CMDLINE_BUILD_TESTS OFF CACHE BOOL "Build the CmdLine unit tests" FORCE)
add_subdirectory(CmdLine)

find_package(GSL REQUIRED)

set(ANALYSISTOOLS_GSL_TARGETS "")
if(TARGET GSL::GSL)
  list(APPEND ANALYSISTOOLS_GSL_TARGETS GSL::GSL)
else()
  if(TARGET GSL::gsl)
    list(APPEND ANALYSISTOOLS_GSL_TARGETS GSL::gsl)
  endif()
  if(TARGET GSL::gslcblas)
    list(APPEND ANALYSISTOOLS_GSL_TARGETS GSL::gslcblas)
  endif()
endif()

if(NOT ANALYSISTOOLS_GSL_TARGETS)
  message(FATAL_ERROR "Could not determine the imported CMake targets provided by GSL.")
endif()

add_library(AnalysisTools
  AnalysisBase.cc
)
add_library(AnalysisTools::AnalysisTools ALIAS AnalysisTools)
set_target_properties(AnalysisTools PROPERTIES
  VERSION ${PROJECT_VERSION}
  SOVERSION ${PROJECT_VERSION_MAJOR}
)

target_compile_features(AnalysisTools PUBLIC cxx_std_20)
target_include_directories(AnalysisTools PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
)
target_link_libraries(AnalysisTools PUBLIC
  CmdLine::CmdLine
  ${ANALYSISTOOLS_GSL_TARGETS}
)

if(ANALYSISTOOLS_BUILD_TESTS)
  add_subdirectory(tests)
endif()

if(ANALYSISTOOLS_BUILD_UNIT_TESTS)
  add_subdirectory(unit-tests)
endif()
