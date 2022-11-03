macro(find_dependencies)

  set(BLA_VENDOR OpenBLAS)

  find_package(BLAS REQUIRED)

  if (NOT BLAS_FOUND)
    message(FATAL_ERROR "OpenBLAS is required with -DWITH_POPSTRAT=ON")
  endif()

  find_package(GSL REQUIRED)
  if (NOT GSL_FOUND)
    message(FATAL_ERROR "GSL is required with -DWITH_POPSTRAT=ON")
  endif()

  set(LAPACK_FOUND TRUE)
  set(LAPACKE_INCDIR "/usr/include/;/usr/include/openblas")
  set(LAPACKE_LIBDIR "/usr/lib64/;/usr/local/lib64")

  if (NOT APPLE)

    find_package(LAPACKE REQUIRED)
    if (NOT LAPACKE_FOUND)
      message(FATAL_ERROR "LAPACKE is requred with -DWITH_POPSTRAT=ON")
    endif()

    if (CONDA_BUILD)
      set(BLAS_INCLUDE_DIRS "BLAS_INCLUDE_DIR-NOTFOUND")
    else()
      find_path(BLAS_INCLUDE_DIRS openblas_config-x86_64.h
        /usr/include
        /usr/include/openblas
        /usr/local/include/
        /usr/local/include/openblas
      )
    endif()

  endif()

  if (BLAS_INCLUDE_DIRS STREQUAL "BLAS_INCLUDE_DIRS-NOTFOUND" OR APPLE)
    set(BLAS_INCLUDE_DIRS "")
    set(OPENBLAS_VERSION "unknown")
  else()
    execute_process(
      COMMAND bash "-c" "grep OPENBLAS_VERSION ${BLAS_INCLUDE_DIRS}/openblas_config-x86_64.h | cut -d' ' -f5"
      OUTPUT_VARIABLE OPENBLAS_VERSION
      OUTPUT_STRIP_TRAILING_WHITESPACE
    )
  endif()

endmacro()

