cmake_minimum_required(VERSION 3.13.0)
project(kmdiff VERSION 1.0.0 LANGUAGES C CXX)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

option(WITH_KMTRICKS "Build kmtricks" ON)
option(WITH_TESTS "Build tests" OFF)
option(WITH_POPSTRAT "Build with population stratification support" ON)
option(WITH_PLUGIN "Build plugins" OFF)

option(CONDA_BUILD "Build inside conda env" OFF)
option(DEV_BUILD "Dev build" OFF)

option(STATIC_BUILD "Static build" OFF)
option(NATIVE_BUILD "Build with -march=native" ON)

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

if (WITH_POPSTRAT)
  include(Depends)
  find_dependencies()
endif()

message(STATUS "CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}")

set(PROJECT_DESCRIPTION "kmdiff - Differential k-mers analysis.")

set(CMAKE_CXX_STANDARD 17)

add_library(build_type_flags INTERFACE)
target_compile_options(build_type_flags INTERFACE
  $<$<CONFIG:Debug>:-O0 -g>
  $<$<CONFIG:Release>:-O3 -DNDEBUG>
  $<$<CONFIG:Profile>:-O3 -ggdb3 -DNDEBUG -fno-inline>
  $<$<CONFIG:Coverage>:-O0 -g --coverage>
)

if (NATIVE_BUILD)
  target_compile_options(build_types_flags INTERFACE "-march=native")
endif()

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
add_library(tests INTERFACE)

if (NOT APPLE AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9.0.0)
  link_libraries(stdc++fs)
endif()

if (NOT MAX_C)
  set(MAX_C 4294967295)
endif()

if (NOT KMER_LIST)
  set (KMER_LIST "32 64 96 128")
endif()

string(REPLACE " " ";" KMER_LIST_ ${KMER_LIST})
string(REPLACE " " "," KMER_LIST ${KMER_LIST})
list(LENGTH KMER_LIST_ KMER_N)

add_compile_definitions(DMAX_C=${MAX_C})

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

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

add_subdirectory(include)
add_subdirectory(src)

if (WITH_PLUGIN)
  add_subdirectory(plugins)
endif()

if (WITH_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()

if (WITH_KMTRICKS)
  add_custom_target(
    copy_km ALL
    DEPENDS KMTRICKS
    COMMAND cp ./thirdparty/kmtricks/bin/kmtricks ${CMAKE_BINARY_DIR}/bin
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  )
  add_dependencies(copy_km KMTRICKS)
endif()

if (WITH_POPSTRAT)
  if (NOT EXISTS ${CMAKE_BINARY_DIR}/bin/smartpca)
    add_custom_target(
      copy_hawk ALL
      DEPENDS HAWK
      COMMAND cp ./thirdparty/HAWK/src/HAWK-build/smartpca ${CMAKE_BINARY_DIR}/bin
      COMMAND cp ../thirdparty/hawk/EIG6.0.1-Hawk/bin/evec2pca.perl ${CMAKE_BINARY_DIR}/bin
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
    add_dependencies(copy_hawk HAWK)
  endif()
  add_custom_target(
    clean_hawk ALL
    DEPENDS HAWK
    COMMAND make clobber
    COMMAND make clean
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/thirdparty/hawk/EIG6.0.1-Hawk/src)
  add_dependencies(clean_hawk HAWK)
endif()

include(GNUInstallDirs)
install (TARGETS ${PROJECT_NAME}-bin DESTINATION ${CMAKE_INSTALL_BINDIR})
install (FILES ${CMAKE_BINARY_DIR}/bin/smartpca
         DESTINATION ${CMAKE_INSTALL_BINDIR}
         PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ)
install (FILES ${CMAKE_BINARY_DIR}/bin/evec2pca.perl
         DESTINATION ${CMAKE_INSTALL_BINDIR}
         PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ)

include(CPackConfig)
