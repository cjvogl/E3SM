#!/bin/bash

# Path to HOMME
export HOMME_ROOT=$HOME/workspace/acme/components/homme

# machine specific settings (compilers, libraries, etc.)
if [[ $HOSTNAME == "cab"* ]]; then
  SYSTEM=cab
  source /usr/local/tools/dotkit/init.sh
  use cmake-3.4.1
  use icc-17.0.174
  use mvapich2-intel-2.2
  use hdf5-intel-parallel-mvapich2-1.10.0
  use netcdf-intel-4.3.3.1
  use netcdf-fortran-intel-4.4.2
elif [[ $HOSTNAME == "quartz"* ]]; then
  SYSTEM=quartz
  module purge
  module load cmake
  module load intel/16.0.3
  module load mvapich2/2.2
  module load hdf5-parallel/1.8.17
  export NETCDF=$HOME/local/netcdf-c-4.3.3.1_quartz_intel_opt
  export HDF5=/usr/tce/packages/hdf5/hdf5-parallel-1.8.17-intel16.0.3-mvapich22.2
fi

# Path to SUNDIALS
export SUNDIALS=$HOME/local/sundials_${SYSTEM}_intel_opt

# HOMME Settings
if [[ $# -gt 1 ]]; then
  NP=$2
  NLEVELS=$3
  NTRACERS=$4
else
  NP=4
  NLEVELS=26
  NTRACERS=4
fi
  
# Set build type ("new" or "update" with default being "update")
if [[ $# -gt 4 ]]; then
  BUILD=$5
else
  BUILD="update"
fi

# Delete and create build directory if necessary
if [[ $BUILD == "new"* ]]; then
  rm -rf $HOMME_ROOT/build_${SYSTEM}
fi
mkdir -p $HOMME_ROOT/build_${SYSTEM} || exit -1
cd $HOMME_ROOT/build_${SYSTEM}

# configure build (if necessary)
if [[ $BUILD == "new" ]]; then
  cmake \
    -C $HOMME_ROOT/cmake/machineFiles/lc.cmake \
    \
    -D PREQX_NP=$NP                  \
    -D PREQX_PLEV=$NLEVELS           \
    -D QSIZE_D=$NTRACERS             \
    \
    $HOMME_ROOT
fi

# build HOMME
make -j$1 theta || exit -1
