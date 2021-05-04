if (NOT BLA_VENDOR)
  set(BLA_VENDOR Generic)
endif()

find_package(BLAS)
find_package(GSL)
find_package(LAPACK)

if (NOT BLAS_FOUND)
  message(FATAL_ERROR "OpenBLAS is required with -DWITH_POPSTRAT=ON")
endif()

if (NOT GSL_FOUND)
  message(FATAL_ERROR "GSL is required with -DWITH_POPSTRAT=ON")
endif()

if (NOT LAPACK_FOUND)
  message(FATAL_ERROR "Lapack is required with -DWITH_POPSTRAT=ON")
endif()