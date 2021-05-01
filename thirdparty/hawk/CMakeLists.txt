cmake_minimum_required(VERSION 3.13.0)
project(hawk)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

include(ExternalProject)
ExternalProject_Add(EIG
    PREFIX EIG
    SOURCE_DIR ${PROJECT_SOURCE_DIR}/EIG6.0.1-Hawk/src
    BUILD_IN_SOURCE 1
    CONFIGURE_COMMAND ""
    BUILD_COMMAND make
    INSTALL_COMMAND "")
