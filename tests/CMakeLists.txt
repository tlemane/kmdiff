file(GLOB TEST_FILES "*_test.cpp")
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/tests/)
SET(CMAKE_CXX_FLAGS "-g -O0 --coverage ${CMAKE_CXX_FLAGS}")
add_executable(${PROJECT_NAME}-tests "main_test.cpp;${PROJECT_SOURCE_DIR}/src/model.cpp;${PROJECT_SOURCE_DIR}/src/simulator.cpp;${PROJECT_SOURCE_DIR}/src/utils.cpp;${PROJECT_SOURCE_DIR}/src/sv.cpp;${PROJECT_SOURCE_DIR}/src/kmtricks_utils.cpp;${PROJECT_SOURCE_DIR}/src/threadpool.cpp" ${TEST_FILES})
add_dependencies(${PROJECT_NAME}-tests ${DEPS})
target_include_directories(${PROJECT_NAME}-tests PUBLIC ${INCLUDES})
target_link_directories(${PROJECT_NAME}-tests PUBLIC ${GTEST_LIB})
target_link_libraries(${PROJECT_NAME}-tests gtest gtest_main pthread)

add_test(
    NAME kmdiff-tests
    COMMAND sh -c "cd ${PROJECT_SOURCE_DIR}/tests/ ; ./${PROJECT_NAME}-tests --verbose"
)

add_custom_target(
    gcov
    COMMAND lcov --directory . --output rapport.info -c
    COMMAND lcov --remove rapport.info -o rapport_filter.info '/usr/**' '${PROJECT_SOURCE_DIR}/thirdparty/**' '${PROJECT_SOURCE_DIR}/tests/**' '${PROJECT_SOURCE_DIR}/build/**' '${PROJECT_SOURCE_DIR}/src/simulator.cpp'
    COMMAND genhtml -o ./rapport_html/ rapport_filter.info 
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/tests/CMakeFiles/kmdiff-tests.dir)
add_dependencies(gcov kmdiff-tests)
