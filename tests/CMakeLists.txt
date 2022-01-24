set(TEST_FILES
    "accumulator_test.cpp"
    "linear_test.cpp"
    "corrector_test.cpp"
    "kmer_test.cpp"
    "kff_test.cpp"
    "factorial_test.cpp"
    "model_test.cpp"
    "utils_test.cpp")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/tests)
add_executable(${PROJECT_NAME}-tests ${TEST_FILES})
target_link_libraries(${PROJECT_NAME}-tests ${PROJECT_NAME}lib headers gtest gtest_main pthread)

add_test(
  NAME kmdiff-tests
  COMMAND sh -c "cd ${PROJECT_SOURCE_DIR}/tests/ ; ./${PROJECT_NAME}-tests --verbose"
)

