cmake_minimum_required(VERSION 3.13.0)
project(alglib VERSION 3.17)

file(GLOB CPP_FILES_ALG RELATIVE ${PROJECT_SOURCE_DIR} "src/*.cpp")

add_library(alglib ${CPP_FILES_ALG})
set_target_properties(alglib PROPERTIES
        ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}"
        LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}"
        PUBLIC_HEADER "src/alglibinternal.h"
)
