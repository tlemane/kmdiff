file(GLOB TEST_FILES "*_test.cpp")
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/tests/)
SET(CMAKE_CXX_FLAGS "-g -O0 ${CMAKE_CXX_FLAGS}")
add_executable(${PROJECT_NAME}-tests "main_test.cpp" ${TEST_FILES})
add_dependencies(${PROJECT_NAME}-tests ${DEPS} kmdifflib)
target_include_directories(${PROJECT_NAME}-tests PUBLIC ${INCLUDES})
target_link_directories(${PROJECT_NAME}-tests PUBLIC ${GTEST_LIB})
target_link_libraries(${PROJECT_NAME}-tests gtest gtest_main pthread kmdifflib)

add_test(
    NAME kmdiff-tests
    COMMAND sh -c "cd ${PROJECT_SOURCE_DIR}/tests/ ; ./${PROJECT_NAME}-tests --verbose"
)
