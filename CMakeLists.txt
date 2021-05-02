cmake_minimum_required(VERSION 3.13.0)
project(kmdiff VERSION 0.1 LANGUAGES C CXX)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

option(DEV_MODE "Dev mode" OFF)
option(COMPILE_TESTS "Compile tests" OFF)
option(BUILD_KMTRICKS_FULL "Build kmtricks binaries" OFF)
option(MAKE_PACKAGE "Build package" OFF)
option(CONDA_BUILD "Build inside conda env" OFF)
option(STATIC "Static build (requires static zlib)" OFF)
option(WITH_POPSIM "Build kmdiff popsim" OFF)
option(WITH_POPSTRAT "Build with population stratification support" ON)

if (WITH_POPSTRAT)
  include(Depends)
endif()

if (STATIC AND NOT APPLE)
  set(CMAKE_EXE_LINKER_FLAGS "-static")
endif()

if (NOT MAX_K)
  set(MAX_K 32)
endif()

if (NOT MAX_C)
  set(MAX_C 4294967295)
endif()

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

message(STATUS "CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}")

set(IS_CONDA_BUILD No)
if (CONDA_BUILD)
  set(IS_CONDA_BUILD Yes)
endif()

set(PROJECT_DESCRIPTION "kmdiff - Differential k-mers analysis.")
set(CONTACT "teo.lemane@inria.fr")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O0 -g -DNDEBUG")

set(DEF_MAX_COUNT 65536)

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
get_sha1(${THIRDPARTY_DIR}/wgsim WGSIM_SHA1)
get_sha1(${THIRDPARTY_DIR}/zlib ZLIB_SHA1)

set(HTSLIB_SHA1 "8cab518f97bbdf43165320c3a4bdef1e9a42893a")
configure_file("${PROJECT_SOURCE_DIR}/include/kmdiff/config.h.in" ${PROJECT_BINARY_DIR}/include/kmdiff/config.hpp)

set(THIRDPARTY_BUILD_DIR ${PROJECT_BINARY_DIR}/thirdparty)

set(SPD_INCLUDE ${THIRDPARTY_DIR}/spdlog/include)
set(ROBIN_INCLUDE ${THIRDPARTY_DIR}/robin-hood-hashing/src/include)
set(KMTRICKS_INCLUDE ${THIRDPARTY_DIR}/kmtricks/libs/ ${THIRDPARTY_DIR}/kmtricks/thirdparty/lz4/lib)
set(FMT_INCLUDE ${THIRDPARTY_DIR}/fmt/include)
set(BCLI_INCLUDE ${THIRDPARTY_DIR}/bcli/include)
set(XXHASH_INCLUDE ${THIRDPARTY_DIR}/xxHash/)
set(GTEST_INCLUDE ${THIRDPARTY_DIR}/googletest/googletest/include)
set(LZ4_INCLUDE ${THIRDPARTY_DIR}/lz4)
set(ALGLIB_INCLUDE ${THIRDPARTY_DIR}/alglib/src)
set(CONFIG_INCLUDE ${PROJECT_BINARY_DIR}/include)
set(KSEQPP_INCLUDE ${THIRDPARTY_DIR}/kseqpp/include ${THIRDPARTY_DIR}/kseqpp/include/kseq++)
set(WGSIM_INCLUDE ${THIRDPARTY_DIR}/wgsim)
set(KFF_INCLUDE ${THIRDPARTY_BUILD_DIR}/KFF/src/KFF-build)
set(HTSLIB_INCLUDE ${THIRDPARTY_BUILD_DIR}/HTSLIB/src/HTSLIB)
set(ZLIB_INCLUDE ${THIRDPARTY_DIR}/zlib)
set(KMTRICKS_LIB ${THIRDPARTY_DIR}/kmtricks/bin/lib)
set(FMT_LIB ${THIRDPARTY_BUILD_DIR}/FMT/src/FMT-build/)
set(XXHASH_LIB ${THIRDPARTY_BUILD_DIR}/XXHASH/src/XXHASH-build)
set(GTEST_LIB ${THIRDPARTY_BUILD_DIR}/GOOGLE_TEST/src/GOOGLE_TEST-build/lib)
set(ALGLIB_LIB ${THIRDPARTY_BUILD_DIR}/ALGLIB/src/ALGLIB-build)
set(WGSIM_LIB ${THIRDPARTY_BUILD_DIR}/WGSIM/src/WGSIM-build)
set(KFF_LIB ${THIRDPARTY_BUILD_DIR}/KFF/src/KFF-build)
set(HTSLIB_LIB ${THIRDPARTY_BUILD_DIR}/HTSLIB/src/HTSLIB)
set(ZLIB_LIB ${THIRDPARTY_BUILD_DIR}/ZLIB/src/ZLIB-build/)

set(INCLUDES ${SPD_INCLUDE} ${KMTRICKS_INCLUDE} ${FMT_INCLUDE} ${BCLI_INCLUDE} ${CONFIG_INCLUDE} ${PROJECT_SOURCE_DIR}/include ${PROJECT_BINARY_DIR}/include ${GTEST_INCLUDE} ${ROBIN_INCLUDE} ${XXHASH_INCLUDE} ${KSEQPP_INCLUDE} ${ALGLIB_INCLUDE} ${WGSIM_INCLUDE} ${KFF_INCLUDE} ${HTSLIB_INCLUDE} ${ZLIB_INCLUDE})

set(DEPS XXHASH SPDLOG FMT KMTRICKS ALGLIB LZ4 KSEQPP WGSIM KFF HTSLIB ZLIB)

include_directories(${INCLUDES})
link_directories(${FMT_LIB} ${XXHASH_LIB} ${KMTRICKS_LIB} ${ALGLIB_LIB} ${WGSIM_LIB} ${KFF_LIB} ${HTSLIB_LIB} ${ZLIB_LIB})

set(HTS_LINK "${HTSLIB_LIB}/libhts.a")
set(ZLIB_LINK "${ZLIB_LIB}/libz.a")

link_libraries(fmt xxhash kmtricks lz4 pthread dl alglib wgsim kff ${HTS_LINK} ${ZLIB_LINK})

add_subdirectory(include)
add_subdirectory(src)

if (COMPILE_TESTS)
  message(STATUS "COMPILE_TESTS=ON - Add target test.")
  enable_testing()
  add_subdirectory(tests)
endif()

if (CMAKE_BUILD_TYPE STREQUAL "Release" AND MAKE_PACKAGE)
  message(STATUS "MAKE_PACKAGE=ON - Add target package.")
  include(CPackConfig)
endif()

add_custom_target(end ALL
        ${CMAKE_COMMAND} -E cmake_echo_color --cyan
        "Binaries: ${PROJECT_SOURCE_DIR}/bin"
        COMMENT "" )

if (COMPILE_TESTS)
  add_dependencies(end ${PROJECT_NAME} ${PROJECT_NAME}-tests)
else()
  add_dependencies(end ${PROJECT_NAME})
endif()

if (BUILD_KMTRICKS_FULL)
  add_custom_target(
    copy_km ALL
    DEPENDS KMTRICKS
    COMMAND cp ./thirdparty/kmtricks/bin/km* ./bin
    COMMAND cp ./thirdparty/kmtricks/kmtricks.py ./bin
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
  add_dependencies(copy_km KMTRICKS)
endif()

if (WITH_POPSTRAT)
  if (NOT EXISTS ${PROJECT_SOURCE_DIR}/bin/smartpca)
    add_custom_target(
      copy_hawk ALL
      DEPENDS HAWK
      COMMAND cp ./thirdparty/hawk/EIG6.0.1-Hawk/src/smartpca ./bin
      COMMAND cp ./thirdparty/hawk/EIG6.0.1-Hawk/bin/evec2pca.perl ./bin
      WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
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
