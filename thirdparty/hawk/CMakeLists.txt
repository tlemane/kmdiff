cmake_minimum_required(VERSION 3.13.0)
project(hawk LANGUAGES C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG")

set(ND "nicksrc")
set(ED "eigensrc")
set(KD "ksrc")
set(SRC "${PROJECT_SOURCE_DIR}/EIG6.0.1-Hawk/src")
set(SRCE "${SRC}/${ED}")
set(SRCK "${SRC}/${KD}")
set(SRCN "${SRC}/${ND}")

if (NOT APPLE)
  if (NOT BLA_VENDOR)
    set(BLA_VENDOR OpenBLAS)
  endif()

  find_package(BLAS)
  find_package(GSL)

  if (NOT BLAS_FOUND)
    message(FATAL_ERROR "OpenBLAS is required with -DWITH_POPSTRAT=ON")
  endif()

  if (NOT GSL_FOUND)
    message(FATAL_ERROR "GSL is required with -DWITH_POPSTRAT=ON")
  endif()
endif()

include(ExternalProject)
ExternalProject_Add(EIGNICK
    PREFIX EIGNICK
    SOURCE_DIR ${PROJECT_SOURCE_DIR}/EIG6.0.1-Hawk/src/nicksrc
    BUILD_IN_SOURCE 1
    CONFIGURE_COMMAND ""
    BUILD_COMMAND make
    INSTALL_COMMAND "")

if (NOT APPLE)
  set(CMAKE_EXE_LINKER_FLAGS ${BLAS_LINKER_FLAGS})
  add_executable(smartpca ${SRCE}/smartpca.c ${SRCE}/eigsubs.c ${SRCE}/exclude.c ${SRCE}/smartsubs.c ${SRCE}/eigx.c ${SRC}/mcio.c ${SRC}/qpsubs.c ${SRC}/admutils.c ${SRC}/egsubs.c ${SRC}/regsubs.c ${SRC}/gval.c ${SRCK}/kjg_fpca.c ${SRCK}/kjg_gsl.c ${SRC}/twsubs.c)
  target_include_directories(smartpca PRIVATE ${SRC}/../include ${BLAS_LIBRAIRIES})
  target_link_directories(smartpca PRIVATE ${BLAS_LIBRARIES} ${SRCN})
  target_link_libraries(smartpca nick openblas GSL::gsl openblas m pthread lapacke lapack)
  add_dependencies(smartpca EIGNICK)
else()
  set(EXT_LIBS /usr/local/opt/gsl/lib /usr/local/opt/openblas/lib /usr/local/opt/lapack/lib)
  set(EXT_INCLUDES /usr/local/opt/gsl/include /usr/local/opt/openblas/include /usr/local/opt/lapack/include)
  add_executable(smartpca ${SRCE}/smartpca.c ${SRCE}/eigsubs.c ${SRCE}/exclude.c ${SRCE}/smartsubs.c ${SRCE}/eigx.c ${SRC}/mcio.c ${SRC}/qpsubs.c ${SRC}/admutils.c ${SRC}/egsubs.c ${SRC}/regsubs.c ${SRC}/gval.c ${SRCK}/kjg_fpca.c ${SRCK}/kjg_gsl.c ${SRC}/twsubs.c)
  target_include_directories(smartpca PRIVATE ${SRC}/../include ${EXT_INCLUDES})
  target_link_directories(smartpca PRIVATE ${EXT_LIBS} ${SRCN})
  target_link_libraries(smartpca nick openblas gsl openblas m pthread lapacke)
  add_dependencies(smartpca EIGNICK)
endif()
