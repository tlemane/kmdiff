#!/bin/bash

function kmdiff_build ()
{
  OPT="-DCMAKE_BUILD_TYPE=${1} -DKMER_LIST=\"${2}\" -DMAX_C=${3} -DWITH_POPSTRAT=${4} -DWITH_TESTS=${5}"
  echo "Options: ${OPT} j=${6}"

  [[ -d "kmdiff_build" ]] || mkdir kmdiff_build
  cd kmdiff_build

  [[ ${8} == 1 ]] && rm CMakeCache.txt

  if [ "$(uname)" != "Darwin" ]; then
    cmake .. -DCMAKE_BUILD_TYPE=${1} -DKMER_LIST="${2}" -DMAX_C=${3} -DWITH_POPSTRAT=${4} -DWITH_TESTS=${5}
  else
    cmake .. -DCMAKE_PREFIX_PATH="/usr/local/opt/openblas;/usr/local/opt/lapack" -DCMAKE_BUILD_TYPE=${1} -DKMER_LIST="${2}" -DMAX_C=${3} -DWITH_POPSTRAT=${4} -DWITH_TESTS=${5}
  fi

  make -j${6}

  if [[ ${7} == 1 ]]; then
    ctest --verbose
  fi
}

function kmdiff_build_conda ()
{
  if [ ! -d "./kmdiff_conda" ]; then
    conda create -y -p ./kmdiff_conda
  fi

  conda activate ./kmdiff_conda

  if [ "$(uname)" == "Darwin" ]; then
    conda install -y -c conda-forge clangxx_osx-64=11.1.0 \
                                    cmake \
                                    zlib \
                                    bzip2 \
                                    openblas \
                                    openblas-devel \
                                    liblapacke \
                                    liblapack \
                                    gsl

    export CC=$(realpath ./kmdiff_conda/bin/x86_64-apple-darwin13.4.0-clang)
    export CXX=$(realpath ./kmdiff_conda/bin/x86_64-apple-darwin13.4.0-clang++)
  else
    conda install -y -c conda-forge gxx_linux-64=9.3.0 \
                                    cmake \
                                    zlib \
                                    bzip2 \
                                    openblas \
                                    openblas-devel \
                                    liblapacke \
                                    liblapack \
                                    gsl

    export CC=$(realpath ./kmdiff_conda/bin/x86_64-conda_cos6-linux-gnu-gcc)
    export CXX=$(realpath ./kmdiff_conda/bin/x86_64-conda_cos6-linux-gnu-g++)
  fi

  OPT="-DCMAKE_BUILD_TYPE=${1} -DKMER_LIST=\"${2}\" -DMAX_C=${3} -DWITH_POPSTRAT=${4}"
  echo "Options: ${OPT} ${6}"

  [[ -d "kmdiff_build" ]] || mkdir kmdiff_build
  cd kmdiff_build

  cmake .. -DCMAKE_PREFIX_PATH=./kmdiff_conda -DKMER_LIST="${2}" -DMAX_C=${3} -DWITH_POPSTRAT=${4} -DWITH_TESTS=${5} -DCMAKE_PREFIX_PATH=$(realpath ../kmdiff_conda)

  make -j${6}

  if [[ ${7} == 1 ]]; then
    ctest --verbose
  fi
}

function usage ()
{
  echo "kmdiff build script."
  echo "Usage: "
  echo "  ./install.sh [-r str] [-k LIST[int]] [-t int] [-c int] [-j int] [-s int] [-p] [-e] [-d] [-h]"
  echo "Options: "
  echo "  -r <Release|Debug> -> build type {Release}."
  echo "  -k <LIST[INT]>     -> k-mer size {\"32 64 96 128\"}."
  echo "  -t <0|1|2>         -> tests: 0 = disabled, 1 = compile, 2 = compile and run {2}."
  echo "  -c <1|2|4>         -> byte per count {4}."
  echo "  -j <INT>           -> nb threads {8}."
  echo "  -s <0|1>           -> population stratification correction 0 = disabled, 1 = enabled {1}"
  echo "                        (-s 1 requires GSL + lapacke + OpenBLAS)"
  echo "  -p                 -> compile with plugins support {disabled}"
  echo "  -e                 -> use conda to install compilers/dependencies {disabled}"
  echo "  -d                 -> delete cmake cache {disabled}"
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
conda=0
cmake_cache=0

while getopts "r:k:t:c:j:s:epdh" option; do
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
      [[ ${tests} == 0 ]] && tests_run=0
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
    e)
      conda=1
      ;;
    d)
      cmake_cache=1
      ;;
    *)
      usage
      ;;
  esac
done

count=$((2**(${count}*8)-1))

if [[ ${conda} -eq 1 ]]; then
  conda_install_path=$(conda info | grep -i 'base environment')
  conda_install_path=$(echo ${conda_install_path} | cut -d' ' -f4)
  source ${conda_install_path}/etc/profile.d/conda.sh
  kmdiff_build_conda ${mode} "${ks}" ${count} ${pop} ${tests_str} ${jopt} ${tests_run} ${cmake_cache}
else
  kmdiff_build ${mode} "${ks}" ${count} ${pop} ${tests_str} ${jopt} ${tests_run} ${cmake_cache}
fi

