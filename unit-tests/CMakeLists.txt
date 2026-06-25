add_executable(unit-tests
  catch_amalgamated.cc
  test_Analysis.cc
  test_GSLRandom.cc
  test_SimpleHist.cc
  test_SimpleHistWithError.cc
)
target_link_libraries(unit-tests PRIVATE AnalysisTools::AnalysisTools)

if(BUILD_TESTING)
  add_test(NAME unit-tests COMMAND unit-tests)
endif()
