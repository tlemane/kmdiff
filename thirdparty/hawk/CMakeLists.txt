cmake_minimum_required(VERSION 3.13.0)
project(hawk LANGUAGES C)

set(CMAKE_C_STANDARD 11)

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

include(Depends)
find_dependencies()

add_library(build_type_flags INTERFACE)
target_compile_options(build_type_flags INTERFACE
  $<$<CONFIG:Debug>:-O0 -g>
  $<$<CONFIG:Release>:-O3 -march=native -DNDEBUG>
)

set(NICK "nicksrc")
set(EIGN "eigensrc")
set(KD "ksrc")

set(SRC_DIR "${CMAKE_SOURCE_DIR}/EIG6.0.1-Hawk/src")

set(SRCN "${SRC_DIR}/${NICK}")
set(SRCE "${SRC_DIR}/${EIGN}")
set(SRCK "${SRC_DIR}/${KD}")

set(SRC_FILES
    ${SRCE}/smartpca.c
    ${SRCE}/eigsubs.c
    ${SRCE}/exclude.c
    ${SRCE}/smartsubs.c
    ${SRCE}/eigx.c
    ${SRCK}/kjg_fpca.c
    ${SRCK}/kjg_gsl.c
    ${SRC_DIR}/mcio.c
    ${SRC_DIR}/qpsubs.c
    ${SRC_DIR}/admutils.c
    ${SRC_DIR}/egsubs.c
    ${SRC_DIR}/regsubs.c
    ${SRC_DIR}/gval.c
    ${SRC_DIR}/twsubs.c
)

include(ExternalProject)
ExternalProject_Add(EIGNICK
  PREFIX EIGNICK
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/EIG6.0.1-Hawk/src/nicksrc
  BUILD_IN_SOURCE 1
  CONFIGURE_COMMAND ""
  BUILD_COMMAND make
  INSTALL_COMMAND ""
)

add_executable(smartpca ${SRC_FILES})
target_include_directories(smartpca PRIVATE ${SRC_DIR}/../include ${BLAS_INCLUDE_DIRS})
target_link_directories(smartpca PRIVATE ${SRCN})
target_link_libraries(smartpca
  nick
  GSL::gsl
  ${BLAS_LIBRARIES}
  m
  pthread
)

if (NOT APPLE)
  target_link_libraries(smartpca ${LAPACKE_LIBRARIES})
endif()
target_compile_options(smartpca PRIVATE "-Wno-alloc-size-larger-than")
add_dependencies(smartpca EIGNICK)
