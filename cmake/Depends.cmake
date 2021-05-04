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