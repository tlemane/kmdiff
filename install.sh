#!/bin/bash

check_conda(){
    which conda &> /dev/null
    if [ $? -ne 0 ]
    then
        echo "$1$reset is required with -e <env>"
        exit 1
    fi
}

function kmdiff_build ()
{
  echo "Options:"
  echo "  -DCMAKE_BUILD_TYPE=${1} -DDEV_MODE=${2} -DWITH_POPSIM=${2} -DBUILD_KMTRICKS_FULL=ON"
  echo -e "  -DMAX_K=${3} -DCOMPILE_TESTS=${4} -DSTATIC=${5} -DMAX_C=${8} -DWITH_POPSTRAT=${9}\n"

  mkdir build
  cd build

  cmake .. -DCMAKE_BUILD_TYPE=${1} \
           -DDEV_MODE=${2} -DWITH_POPSIM=${2} \
           -DBUILD_KMTRICKS_FULL=ON \
           -DMAX_K=${3} \
           -DCOMPILE_TESTS=${4} \
           -DSTATIC=${5} \
           -DMAX_C=${8} \
           -DWITH_POPSTRAT=${9}
  make -j${6}

  if [[ ${7} == 1 ]]; then
    ctest --verbose
  fi
}

function conda_build ()
{
  conda create -y -p ${10}
  conda activate ${10}
  conda install -y mamba -c conda-forge
  mamba install -y liblapack liblapacke gsl openblas zlib bzip2 -c conda-forge

  if [ "$(uname)" == "Darwin" ]; then
    mamba install -y clangxx_osx-64>=11.1.0 -c conda-forge
  elif [ "$(expr substr $(uname -s) 1 5)" == "Linux" ]; then
    mamba install -y gxx_linux-64>=9.3.0 -c conda-forge
    export CXX=${10}/bin/x86_64-conda_cos6-linux-gnu-c++
    export CC=${10}/bin/x86_64-conda_cos6-linux-gnu-gcc
  fi

  echo "Options:"
  echo "  -DCMAKE_BUILD_TYPE=${1} -DDEV_MODE=${2} -DWITH_POPSIM=${2} -DBUILD_KMTRICKS_FULL=ON"
  echo -e "  -DMAX_K=${3} -DCOMPILE_TESTS=${4} -DSTATIC=${5} -DMAX_C=${8} -DWITH_POPSTRAT=${9}\n"

  mkdir build
  cd build

  cmake .. -DCMAKE_BUILD_TYPE=${1} \
           -DDEV_MODE=${2} -DWITH_POPSIM=${2} \
           -DBUILD_KMTRICKS_FULL=ON \
           -DMAX_K=${3} \
           -DCOMPILE_TESTS=${4} \
           -DSTATIC=${5} \
           -DMAX_C=${8} \
           -DWITH_POPSTRAT=${9} \
           -DCMAKE_PREFIX_PATH=${CONDA_PREFIX}
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
  echo "  -t <0|1|2>         -> tests: 0 = disabled, 1 = compile, 2 = compile and run {2}."
  echo "  -c <1|2|4>         -> byte per count {4}."
  echo "  -j <INT>           -> nb threads {8}."
  echo "  -p <0|1>           -> population stratification correction 0 = disabled, 1 = enabled {1}"
  echo "                        (-p 1 requires GSL + lapack + OpenBLAS)"
  echo "  -e <STR>           -> build inside a conda environment"
  echo "  -d                 -> dev build {disabled}."
  echo "  -s                 -> static build {disabled}."
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
conda_env=""


while getopts "r:k:t:c:j:p:e:dsh" option; do
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
    e)
      conda_env=${OPTARG}
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
      pop=${OPTARG}
      [[ ${pop} == 0 ]] || [[ ${pop} == 1 ]] || usage
      [[ ${pop} == 0 ]] && pop="OFF"
      ;;
    *)
      usage
      ;;
  esac
done

count=$((2**(${count}*8)-1))

if [[ ${conda_env} == "" ]]; then
  kmdiff_build ${mode} ${dev} ${ks} ${tests_str} ${static} ${jopt} ${tests_run} ${count} ${pop}
else
  check_conda
  conda_env=$(realpath ${conda_env})
  conda_install_path=$(conda info | grep -i 'base environment')
  conda_install_path=$(echo ${conda_install_path} | cut -d' ' -f4)
  source ${conda_install_path}/etc/profile.d/conda.sh
  conda_build ${mode} ${dev} ${ks} ${tests_str} ${static} ${jopt} ${tests_run} ${count} ${pop} ${conda_env}
fi