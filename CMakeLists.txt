cmake_minimum_required(VERSION 3.13.0)
project(kmdiff VERSION 0.0.1 LANGUAGES C CXX)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

option(DEV_MODE "Dev mode" OFF)
option(WITH_TESTS "Compile tests" OFF)
option(MAKE_PACKAGE "Build package" OFF)
option(CONDA_BUILD "Build inside conda env" OFF)
option(STATIC "Static build" OFF)
option(WITH_POPSTRAT "Build with population stratification support" ON)
option(WITH_PLUGIN "Build plugins" OFF)

#set(CMAKE_FIND_LIBRARY_SUFFIXES ".a" ".so" ".dylib")

if (WITH_POPSTRAT)
  include(Depends)
  find_dependencies()
endif()

if (STATIC AND NOT APPLE)
  set(CMAKE_EXE_LINKER_FLAGS "-static")
endif()

if (NOT MAX_C)
  set(MAX_C 4294967295)
endif()

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

message(STATUS "CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}")

set(PROJECT_DESCRIPTION "kmdiff - Differential k-mers analysis.")

set(CMAKE_CXX_STANDARD 17)

add_library(build_type_flags INTERFACE)
target_compile_options(build_type_flags INTERFACE
  $<$<CONFIG:Debug>:-O0 -g>
  $<$<CONFIG:Release>:-O3 -march=native -DNDEBUG>
  $<$<CONFIG:Profile>:-O3 -ggdb3 -DNDEBUG -fno-inline>
  $<$<CONFIG:Coverage>:-O0 -g --coverage>
)

add_library(warning_flags INTERFACE)
target_compile_options(warning_flags INTERFACE
  "-Wall"
  "-Wextra"
  "-Wno-char-subscripts"
  "-Wno-unused-parameter"
)

add_library(headers INTERFACE)
target_include_directories(headers INTERFACE
  ${CMAKE_SOURCE_DIR}/include
  ${CMAKE_BINARY_DIR}/include
)
add_library(links INTERFACE)
target_link_libraries(links INTERFACE pthread dl)
add_library(deps INTERFACE)

if (NOT KMER_LIST)
  set (KMER_LIST "32 64 96 128")
endif()

string(REPLACE " " ";" KMER_LIST_ ${KMER_LIST})
string(REPLACE " " "," KMER_LIST ${KMER_LIST})
list(LENGTH KMER_LIST_ KMER_N)

add_compile_definitions(DMAX_C=${MAX_C})
if (DEV_MODE)
  add_compile_definitions(KMDIFF_DEV_MODE)
endif()

if (WITH_POPSIM)
  add_compile_definitions(WITH_POPSIM)
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)

add_subdirectory(thirdparty)

set(THIRDPARTY_DIR ${PROJECT_SOURCE_DIR}/thirdparty)

include(GitUtils)
get_sha1(${PROJECT_SOURCE_DIR} GIT_SHA1)
get_sha1(${THIRDPARTY_DIR}/bcli BCLI_SHA1)
get_sha1(${THIRDPARTY_DIR}/kmtricks KMTRICKS_SHA1)
get_sha1(${THIRDPARTY_DIR}/fmt FMT_SHA1)
get_sha1(${THIRDPARTY_DIR}/kff-cpp-api KFF_SHA1)
get_sha1(${THIRDPARTY_DIR}/lz4 LZ4_SHA1)
get_sha1(${THIRDPARTY_DIR}/robin-hood-hashing ROBIN_SHA1)
get_sha1(${THIRDPARTY_DIR}/spdlog SPDLOG_SHA1)
get_sha1(${THIRDPARTY_DIR}/xxHash XXHASH_SHA1)
get_sha1(${THIRDPARTY_DIR}/indicators INDICATORS_SHA1)

configure_file("${PROJECT_SOURCE_DIR}/include/kmdiff/config.h.in" ${PROJECT_BINARY_DIR}/include/kmdiff/config.hpp)

set(THIRDPARTY_BUILD_DIR ${PROJECT_BINARY_DIR}/thirdparty)

add_subdirectory(include)
add_subdirectory(src)

if (WITH_PLUGIN)
  add_subdirectory(plugins)
endif()

#if (WITH_TESTS)
#  enable_testing()
#  add_subdirectory(tests)
#endif()


#if (CMAKE_BUILD_TYPE STREQUAL "Release" AND MAKE_PACKAGE)
#  message(STATUS "MAKE_PACKAGE=ON - Add target package.")
#  include(CPackConfig)
#endif()
#
#add_custom_target(end ALL
#        ${CMAKE_COMMAND} -E cmake_echo_color --cyan
#        "Binaries: ${PROJECT_SOURCE_DIR}/bin"
#        COMMENT "" )
#
#if (COMPILE_TESTS)
#  add_dependencies(end ${PROJECT_NAME} ${PROJECT_NAME}-tests)
#else()
#  add_dependencies(end ${PROJECT_NAME})
#endif()
#

add_custom_target(
  copy_km ALL
  DEPENDS KMTRICKS
  COMMAND cp ./thirdparty/kmtricks/bin/kmtricks ./bin
  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
)
add_dependencies(copy_km KMTRICKS)

if (WITH_POPSTRAT)
  if (NOT EXISTS ${PROJECT_SOURCE_DIR}/bin/smartpca)
    add_custom_target(
      copy_hawk ALL
      DEPENDS HAWK
      COMMAND cp ./thirdparty/HAWK/src/HAWK-build/smartpca ../bin
      COMMAND cp ../thirdparty/hawk/EIG6.0.1-Hawk/bin/evec2pca.perl ../bin
      WORKING_DIRECTORY ${PROJECT_BINARY_DIR})
    add_dependencies(copy_hawk HAWK)
  endif()
  add_custom_target(
    clean_hawk ALL
    DEPENDS HAWK
    COMMAND make clobber
    COMMAND make clean
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/thirdparty/hawk/EIG6.0.1-Hawk/src)
  add_dependencies(clean_hawk HAWK)
endif()
