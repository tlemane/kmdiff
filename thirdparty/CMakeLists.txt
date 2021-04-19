include(ExternalProject)

ExternalProject_Add(FMT
        PREFIX FMT
        SOURCE_DIR ${PROJECT_SOURCE_DIR}/thirdparty/fmt
        INSTALL_COMMAND "")

ExternalProject_Add(KMTRICKS
        PREFIX KMTRICKS
        SOURCE_DIR ${PROJECT_SOURCE_DIR}/thirdparty/kmtricks
        INSTALL_COMMAND "")

ExternalProject_Add(GOOGLE_TEST
        PREFIX GOOGLE_TEST
        SOURCE_DIR ${PROJECT_SOURCE_DIR}/thirdparty/googletest
        INSTALL_COMMAND "")

ExternalProject_Add(SPDLOG
        PREFIX SPDLOG
        SOURCE_DIR ${PROJECT_SOURCE_DIR}/thirdparty/spdlog
        INSTALL_COMMAND "")

ExternalProject_Add(XXHASH
        PREFIX XXHASH
        SOURCE_DIR ${PROJECT_SOURCE_DIR}/thirdparty/xxHash
        SOURCE_SUBDIR "cmake_unofficial"
        CMAKE_ARGS -DBUILD_SHARED_LIBS=OFF
        INSTALL_COMMAND "")

ExternalProject_Add(LZ4
        PREFIX LZ4
        SOURCE_DIR ${PROJECT_SOURCE_DIR}/thirdparty/lz4
        SOURCE_SUBDIR build/cmake
        CMAKE_ARGS -DBUILD_STATIC_LIBS=ON -DBUILD_SHARED_LIBS=OFF
        INSTALL_COMMAND "")

ExternalProject_Add(ALGLIB
        PREFIX ALGLIB
        SOURCE_DIR ${PROJECT_SOURCE_DIR}/thirdparty/alglib
        INSTALL_COMMAND "")

ExternalProject_Add(KFF
        PREFIX KFF
        SOURCE_DIR ${PROJECT_SOURCE_DIR}/thirdparty/kff-cpp-api
        INSTALL_COMMAND "")

ExternalProject_Add(KSEQPP
        PREFIX KSEQPP
        SOURCE_DIR ${PROJECT_SOURCE_DIR}/thirdparty/kseqpp
        INSTALL_COMMAND "")

ExternalProject_Add(WGSIM
        PREFIX WGSIM
        SOURCE_DIR ${PROJECT_SOURCE_DIR}/thirdparty/wgsim
        INSTALL_COMMAND "")

ExternalProject_Add(HTSLIB
        PREFIX HTSLIB
        URL https://github.com/samtools/htslib/releases/download/1.12/htslib-1.12.tar.bz2
        BUILD_IN_SOURCE 1
        CONFIGURE_COMMAND ${PROJECT_BINARY_DIR}/thirdparty/HTSLIB/src/HTSLIB/configure --disable-libcurl --disable-bz2 --disable-lzma
        BUILD_COMMAND make
        INSTALL_COMMAND "")