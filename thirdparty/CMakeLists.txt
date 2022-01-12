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
)
target_include_directories(headers INTERFACE ${THIRD_DIR}/fmt/include)
target_link_directories(links INTERFACE ${THIRD_BINDIR}/FMT/src/FMT-build)
target_link_libraries(links INTERFACE fmt)
add_dependencies(deps FMT)

ExternalProject_Add(
  KMTRICKS
  PREFIX KMTRICKS
  SOURCE_DIR ${PROJECT_SOURCE_DIR}/thirdparty/kmtricks
  INSTALL_COMMAND ""
)
target_include_directories(headers INTERFACE ${THIRD_DIR}/kmtricks/libs)
target_link_directories(links INTERFACE ${THIRD_DIR}/kmtricks/bin/lib)
target_link_libraries(links INTERFACE kmtricks)

ExternalProject_Add(GOOGLE_TEST
        PREFIX GOOGLE_TEST
        SOURCE_DIR ${PROJECT_SOURCE_DIR}/thirdparty/googletest
        INSTALL_COMMAND "")

ExternalProject_Add(
  SPDLOG
  PREFIX SPDLOG
  SOURCE_DIR ${PROJECT_SOURCE_DIR}/thirdparty/spdlog
  INSTALL_COMMAND ""
)
target_include_directories(headers INTERFACE ${THIRD_DIR}/spdlog/include)
add_dependencies(deps SPDLOG)

ExternalProject_Add(
  XXHASH
  PREFIX XXHASH
  SOURCE_DIR ${PROJECT_SOURCE_DIR}/thirdparty/xxHash
  SOURCE_SUBDIR "cmake_unofficial"
  CMAKE_ARGS -DBUILD_SHARED_LIBS=OFF
  INSTALL_COMMAND ""
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
)
target_include_directories(headers INTERFACE ${THIRD_BINDIR}/KFF/src/KFF-build)
target_link_directories(links INTERFACE ${THIRD_BINDIR}/KFF/src/KFF-build)
target_link_libraries(links INTERFACE kff)
add_dependencies(deps KFF)

ExternalProject_Add(
  KSEQPP
  PREFIX KSEQPP
  SOURCE_DIR ${PROJECT_SOURCE_DIR}/thirdparty/kseqpp
  INSTALL_COMMAND ""
)
target_include_directories(headers INTERFACE ${THIRD_DIR}/kseqpp/include ${THIRD_DIR}/kseqpp/include/kseq++)
add_dependencies(deps KSEQPP)

ExternalProject_Add(
  WGSIM
  PREFIX WGSIM
  SOURCE_DIR ${PROJECT_SOURCE_DIR}/thirdparty/wgsim
  INSTALL_COMMAND ""
)
target_include_directories(headers INTERFACE ${THIRD_DIR}/wgsim)
target_link_directories(links INTERFACE ${THIRD_BINDIR}/WGSIM/src/WGSIM-build)
target_link_libraries(links INTERFACE wgsim)
add_dependencies(deps WGSIM)


ExternalProject_Add(
  HTSLIB
  PREFIX HTSLIB
  URL https://github.com/samtools/htslib/releases/download/1.12/htslib-1.12.tar.bz2
  BUILD_IN_SOURCE 1
  BUILD_BYPRODUCTS ${PROJECT_BINARY_DIR}/thirdparty/HTSLIB/src/HTSLIB/libhts.a
  CONFIGURE_COMMAND ${PROJECT_BINARY_DIR}/thirdparty/HTSLIB/src/HTSLIB/configure --disable-libcurl --disable-bz2 --disable-lzma
  BUILD_COMMAND make -j4
  INSTALL_COMMAND ""
)
target_include_directories(headers INTERFACE ${THIRD_DIR}/htslib)

if (WITH_POPSTRAT)
        ExternalProject_Add(HAWK
                PREFIX HAWK
                SOURCE_DIR ${PROJECT_SOURCE_DIR}/thirdparty/hawk
                INSTALL_COMMAND "")
endif()

