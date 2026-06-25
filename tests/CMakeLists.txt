add_executable(test-all-types
  test-all-types.cc
)
target_link_libraries(test-all-types PRIVATE AnalysisTools::AnalysisTools)

add_executable(test-map-speed
  test-map-speed.cc
)
target_link_libraries(test-map-speed PRIVATE AnalysisTools::AnalysisTools)
