include_directories(catch2) # v2.2.1

### Unit tests
file(GLOB_RECURSE FILENAMES *.cpp)
list(FILTER FILENAMES EXCLUDE REGEX "/test/algorithms/quality.cpp")

add_executable(run_tests ${FILENAMES})
target_link_libraries(run_tests PUBLIC mockturtle)
if (ENABLE_COVERAGE)
  target_link_libraries(run_tests PUBLIC gcov)
endif()
target_compile_definitions(run_tests PUBLIC CATCH_CONFIG_CONSOLE_WIDTH=300)

### Quality tests
add_executable(quality_tests test.cpp algorithms/quality.cpp)
target_link_libraries(quality_tests PUBLIC mockturtle)
target_compile_definitions(quality_tests PUBLIC BENCHMARKS_PATH="${CMAKE_CURRENT_SOURCE_DIR}/../experiments/benchmarks")
target_compile_definitions(quality_tests PUBLIC CATCH_CONFIG_CONSOLE_WIDTH=300)
