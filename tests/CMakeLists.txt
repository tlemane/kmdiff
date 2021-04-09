file(GLOB TEST_FILES "*_test.cpp")
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/tests/)
add_executable(${PROJECT_NAME}-tests "main_test.cpp" ${TEST_FILES})
add_dependencies(${PROJECT_NAME}-tests ${DEPS} kmsv)
target_include_directories(${PROJECT_NAME}-tests PUBLIC ${INCLUDES})
target_link_directories(${PROJECT_NAME}-tests PUBLIC ${GTEST_LIB})
target_link_libraries(${PROJECT_NAME}-tests gtest gtest_main pthread kmsv)

add_test(
    NAME kmtricks_sv_tests
    COMMAND sh -c "cd ${PROJECT_SOURCE_DIR}/tests/ ; ./${PROJECT_NAME}-tests --verbose"
)
