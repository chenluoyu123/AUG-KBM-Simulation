#!/bin/bash


export jlchen_gacode_third_party_root=/public/home/jianx/mycode/third_party
module purge
#module load  compiler/devtoolset/7.3.1   mpi/hpcx/2.11.0/gcc-7.3.1   compiler/rocm/dtk-22.10.1
module use ${jlchen_gacode_third_party_root}/modulefiles
module load  compiler/devtoolset/7.3.1   openmpi-MINGYUE
module load zlib-jlchen
#module load openmpi-jlchen
module load hdf5-jlchen
module load netcdf-jlchen
module load fftw-jlchen

export GACODE_PLATFORM=MINGYUE # set this for compile
export GACODE_ROOT=/public/home/jianx/mycode/gacode_2021_gfortran
. ${GACODE_ROOT}/shared/bin/gacode_setup


echo ${GACODE_ROOT}


$GACODE_ROOT/shared/bin/gacode_qsub  -code "cgyro" -queue hfacnormal01 -w 12:00:00 -n 1024 -s
