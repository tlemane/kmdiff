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

  find_package(LAPACKE REQUIRED)
  if (NOT LAPACKE_FOUND)
    message(FATAL_ERROR "LAPACKE is requred with -DWITH_POPSTRAT=ON")
  endif()

  find_path(BLAS_INCLUDE_DIRS openblas_config-x86_64.h
    /usr/include
    /usr/include/openblas
    /usr/local/include/
    /usr/local/include/openblas
  )

  execute_process(
    COMMAND bash "-c" "grep OPENBLAS_VERSION ${BLAS_INCLUDE_DIRS}/openblas_config-x86_64.h | cut -d' ' -f5"
    OUTPUT_VARIABLE OPENBLAS_VERSION
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
endmacro()

