#!/bin/bash

function kmdiff_build ()
{
  echo "mode=${1} dev=${2} k=${3} c=${8} tests=${4} static=${5} j=${6} run=${7} pop=${9}"

  mkdir build
  cd build
  cmake .. -DCMAKE_BUILD_TYPE=${1} \
           -DDEV_MODE=${2} -DWITH_POPSIM=${2} \
           -DBUILD_KMTRICKS_FULL=ON \
           -DMAX_K=${3}
           -DCOMPILE_TESTS=${4} \
           -DSTATIC=${5} \
           -DMAX_C=${8} \
           -DWITH_POPSTRAT=${9}
  make -j${6}

  if [[ ${7} == 1 ]]; then
    ctest --verbose
  fi
}

function usage ()
{
  echo "kmdiff build script."
  echo "Usage: "
  echo "  ./install.sh [-r str] [-k int] [-t int] [-c int] [-p int] [-j int] [-d] [-s] [-h]"
  echo "Options: "
  echo "  -r <Release|Debug> -> build type {Release}."
  echo "  -k <0|32|64>       -> max k-mer size, 0 = both {0}."
  echo "  -t <0|1|2>         -> tests: 0 = disable, 1 = compile, 2 = compile and run {2}."
  echo "  -c <1|2|4>         -> byte per count {4}."
  echo "  -j <INT>           -> nb threads {8}."
  echo "  -p <0|1>           -> population stratification correction 0 = disable, 1 = enable {1}"
  echo "                        (-p 1 requires GSL + lapack + OpenBLAS)"
  echo "  -d                 -> dev build {disable}."
  echo "  -s                 -> static build {disable}."
  echo "  -h                 -> show help."
  exit 1
}

rel="Release"
deb="Debug"

mode="Release"
ks="ALL"
count=4
dev="OFF"
static="OFF"
tests=2
tests_str="ON"
tests_run=1
jopt=8
pop="ON"

while getopts "r:k:t:c:j:p:dsh" option; do
  case "$option" in
    r)
      mode=${OPTARG}
      [[ ${mode} != ${rel} && ${mode} != ${deb} ]] && usage
      ;;
    k)
      ks=${OPTARG}
      [[ ${ks} == "32" ]] || [[ ${ks} == "64" ]] || [[ ${ks} == "0" ]] || usage
      [[ ${ks} == "0" ]] && ks="ALL"
      ;;
    t)
      tests=${OPTARG}
      [[ ${tests} == 0 ]] || [[ ${tests} == 1 ]] || [[ ${tests} == 2 ]] || usage
      [[ ${tests} == 0 ]] && tests_str="OFF"
      [[ ${tests} == 1 ]] && tests_run=0
      ;;
    c)
      count=${OPTARG}
      [[ ${count} == 1 ]] || [[ ${count} == 2 ]] || [[ ${count} == 4 ]] || usage
      ;;
    j)
      jopt=${OPTARG}
      ;;
    d)
      dev="ON"
      ;;
    s)
      static="ON"
      ;;
    h)
      usage
      ;;
    p)
      [[ ${pop} == 0 ]] || [[ ${pop} == 1 ]] || usage
      [[ ${pop} == 0 ]] && pop="OFF"
      ;;
    *)
      usage
      ;;
  esac
done

count=$((2**(${count}*8)-1))
echo ${count}
exit 1
kmdiff_build ${mode} ${dev} ${ks} ${tests_str} ${static} ${jopt} ${tests_run} ${count} ${pop}
