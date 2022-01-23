#!/bin/bash

function kmdiff_build ()
{
  OPT="-DCMAKE_BUILD_TYPE=${1} -DKMER_LIST=${2} -DMAX_C=${3} -DWITH_POPSTRAT=${4}"
  echo "Options: ${OPT} ${6}"

  mkdir kmdiff_build
  cd kmdiff_build

  cmake .. ${OPT}

  make -j${6}

  if [[ ${7} == 1 ]]; then
    ctest --verbose
  fi
}

function usage ()
{
  echo "kmdiff build script."
  echo "Usage: "
  echo "  ./install.sh [-r str] [-k LIST[int]] [-t int] [-c int] [-j int] [-s int] [-p] [-h]"
  echo "Options: "
  echo "  -r <Release|Debug> -> build type {Release}."
  echo "  -k <LIST[INT]>     -> k-mer size {\"32 64 96 128\"}."
  echo "  -t <0|1|2>         -> tests: 0 = disabled, 1 = compile, 2 = compile and run {2}."
  echo "  -c <1|2|4>         -> byte per count {4}."
  echo "  -j <INT>           -> nb threads {8}."
  echo "  -s <0|1>           -> population stratification correction 0 = disabled, 1 = enabled {1}"
  echo "                        (-s 1 requires GSL + lapacke + OpenBLAS)"
  echo "  -p                 -> compile with plugins support {disabled}"
  echo "  -h                 -> show help."
  exit 1
}

rel="Release"
deb="Debug"

ks="32 64 96 128"
mode="Release"
count=4
tests=2
tests_str="ON"
tests_run=1
jopt=8
pop="ON"
plugin="OFF"

while getopts "r:k:t:c:j:s:ph" option; do
  case "$option" in
    r)
      mode=${OPTARG}
      [[ ${mode} != ${rel} && ${mode} != ${deb} ]] && usage
      ;;
    k)
      ks=${OPTARG}
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
    h)
      usage
      ;;
    s)
      pop=${OPTARG}
      [[ ${pop} == 0 ]] || [[ ${pop} == 1 ]] || usage
      [[ ${pop} == 0 ]] && pop="OFF"
      ;;
    p)
      plugin="OFF"
      ;;
    *)
      usage
      ;;
  esac
done

count=$((2**(${count}*8)-1))

kmdiff_build ${mode} "${ks}" ${count} ${pop} ${test_str} ${jopt} ${tests_run}
