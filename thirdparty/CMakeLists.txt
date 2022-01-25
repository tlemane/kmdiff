include(ExternalProject)

find_package(ZLIB REQUIRED)
target_link_libraries(links INTERFACE ZLIB::ZLIB)

set(THIRD_DIR ${CMAKE_SOURCE_DIR}/thirdparty)
set(THIRD_BINDIR ${CMAKE_BINARY_DIR}/thirdparty)

target_include_directories(headers INTERFACE
  ${THIRD_DIR}/bcli/include
  ${THIRD_DIR}/robin-hood-hashing/src/include
)

ExternalProject_Add(
  FMT
  PREFIX FMT
  SOURCE_DIR ${PROJECT_SOURCE_DIR}/thirdparty/fmt
  INSTALL_COMMAND ""
  LOG_CONFIGURE ON
  LOG_BUILD ON
)
target_include_directories(headers SYSTEM INTERFACE ${THIRD_DIR}/fmt/include)
target_link_directories(links INTERFACE ${THIRD_BINDIR}/FMT/src/FMT-build)
target_link_libraries(links INTERFACE fmt)
add_dependencies(deps FMT)

string(REPLACE "," " " KMER_LIST_KM ${KMER_LIST})

if (WITH_KMTRICKS)
  ExternalProject_Add(
    KMTRICKS
    PREFIX KMTRICKS
    SOURCE_DIR ${PROJECT_SOURCE_DIR}/thirdparty/kmtricks
    CMAKE_ARGS "-DKMER_LIST=${KMER_LIST_KM}" "-DMAX_C=${MAX_C}"
    INSTALL_COMMAND ""
    LOG_CONFIGURE ON
    LOG_BUILD ON
  )
else()
  ExternalProject_Add(
    KMTRICKS
    PREFIX KMTRICKS
    SOURCE_DIR ${PROJECT_SOURCE_DIR}/thirdparty/kmtricks
    CMAKE_ARGS "-DKMER_LIST=${KMER_LIST_KM}" "-DMAX_C=${MAX_C}"
    INSTALL_COMMAND ""
    BUILD_COMMAND ""
    LOG_CONFIGURE ON
    LOG_BUILD ON
  )
endif()

target_include_directories(headers SYSTEM INTERFACE
  ${THIRD_DIR}/kmtricks/include
  ${THIRD_DIR}/kmtricks/thirdparty/TurboPFor-Integer-Compression
  ${THIRD_BINDIR}/KMTRICKS/src/KMTRICKS-build/include
)

ExternalProject_Add(
  GOOGLE_TEST
  PREFIX GOOGLE_TEST
  SOURCE_DIR ${PROJECT_SOURCE_DIR}/thirdparty/googletest
  INSTALL_COMMAND ""
  LOG_CONFIGURE ON
  LOG_BUILD ON
)
target_include_directories(tests SYSTEM INTERFACE ${THIRD_DIR}/googletest/googletest/include)
target_link_directories(tests INTERFACE ${THIRD_BINDIR}/GOOGLE_TEST/src/GOOGLE_TEST-build/lib)
target_link_libraries(tests INTERFACE gtest gtest_main pthread)
add_dependencies(deps GOOGLE_TEST)

ExternalProject_Add(
  SPDLOG
  PREFIX SPDLOG
  SOURCE_DIR ${PROJECT_SOURCE_DIR}/thirdparty/spdlog
  INSTALL_COMMAND ""
  LOG_CONFIGURE ON
  LOG_BUILD ON
)
target_include_directories(headers SYSTEM INTERFACE ${THIRD_DIR}/spdlog/include)
add_dependencies(deps SPDLOG)

#Set standard to C99 because of bug in xxHash related to C11 static_assert
ExternalProject_Add(
  XXHASH
  PREFIX XXHASH
  SOURCE_DIR ${PROJECT_SOURCE_DIR}/thirdparty/xxHash
  SOURCE_SUBDIR "cmake_unofficial"
  CMAKE_ARGS -DBUILD_SHARED_LIBS=OFF -DXXHASH_BUILD_XXHSUM=0FF -DCMAKE_C_STANDARD=99
  INSTALL_COMMAND ""
  LOG_CONFIGURE ON
  LOG_BUILD ON
)
target_include_directories(headers INTERFACE ${THIRD_DIR}/xxHash)
target_link_directories(links INTERFACE ${THIRD_BINDIR}/XXHASH/src/XXHASH-build)
target_link_libraries(links INTERFACE xxhash)
add_dependencies(deps XXHASH)

ExternalProject_Add(
  LZ4
  PREFIX LZ4
  SOURCE_DIR ${PROJECT_SOURCE_DIR}/thirdparty/lz4
  SOURCE_SUBDIR build/cmake
  CMAKE_ARGS -DBUILD_STATIC_LIBS=ON -DBUILD_SHARED_LIBS=OFF
  INSTALL_COMMAND ""
  LOG_CONFIGURE ON
  LOG_BUILD ON
)
target_include_directories(headers INTERFACE ${THIRD_DIR}/lz4/lib)
target_link_directories(links INTERFACE ${THIRD_BINDIR}/LZ4/src/LZ4-build)
target_link_libraries(links INTERFACE lz4)
add_dependencies(deps LZ4)

ExternalProject_Add(
  ALGLIB
  PREFIX ALGLIB
  SOURCE_DIR ${PROJECT_SOURCE_DIR}/thirdparty/alglib
  INSTALL_COMMAND ""
  LOG_CONFIGURE ON
  LOG_BUILD ON
)
target_include_directories(headers INTERFACE ${THIRD_DIR}/alglib/src)
target_link_directories(links INTERFACE ${THIRD_BINDIR}/ALGLIB/src/ALGLIB-build)
target_link_libraries(links INTERFACE alglib)
add_dependencies(deps ALGLIB)

ExternalProject_Add(
  KFF
  PREFIX KFF
  SOURCE_DIR ${PROJECT_SOURCE_DIR}/thirdparty/kff-cpp-api
  INSTALL_COMMAND ""
  LOG_CONFIGURE ON
  LOG_BUILD ON
)
target_include_directories(headers SYSTEM INTERFACE ${THIRD_BINDIR}/KFF/src/KFF-build)
target_link_directories(links INTERFACE ${THIRD_BINDIR}/KFF/src/KFF-build)
target_link_libraries(links INTERFACE kff)
add_dependencies(deps KFF)

ExternalProject_Add(
  KSEQPP
  PREFIX KSEQPP
  SOURCE_DIR ${PROJECT_SOURCE_DIR}/thirdparty/kseqpp
  INSTALL_COMMAND ""
  LOG_CONFIGURE ON
  LOG_BUILD ON
)
target_include_directories(headers INTERFACE ${THIRD_DIR}/kseqpp/include ${THIRD_DIR}/kseqpp/include/kseq++)
add_dependencies(deps KSEQPP)

target_include_directories(headers INTERFACE ${THIRD_DIR}/indicators/include)

if (WITH_POPSTRAT)
  ExternalProject_Add(HAWK
    PREFIX HAWK
    SOURCE_DIR ${CMAKE_SOURCE_DIR}/thirdparty/hawk
    CMAKE_ARGS -DCMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH}
    INSTALL_COMMAND ""
    LOG_CONFIGURE ON
    LOG_BUILD ON
  )
endif()

