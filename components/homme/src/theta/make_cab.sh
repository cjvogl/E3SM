#!/bin/bash

source /usr/local/tools/dotkit/init.sh

use icc-17.0.174
use mvapich2-intel-2.2
use netcdf-intel-4.3.3.1
use netcdf-fortran-intel-4.4.2

export NP=4
export PLEV=26
export ACME_PATH=~/workspace/acme
export NETCDF_PATH=/usr/local/tools/netcdf-fortran-intel-4.4.2
export HOMME_PATH=${ACME_PATH}/components/homme
export CONFIGH_PATH=${HOMME_PATH}/build/src/theta

export PIO_PATH=${HOMME_PATH}/utils/cime/externals/pio1/pio
export PIO_PATH2=${HOMME_PATH}/build/utils/cime/externals/pio1/pio
export GPTL_PATH=${HOMME_PATH}/build/utils/cime/share/timing
export SHR_PATH=${HOMME_PATH}/utils/csm_share

export FC=mpif90
export DEFINES='-DCPRINTEL -DHAVE_CONFIG_H -DINCLUDE_CMAKE_FCI -DSPMD'
export FFLAGS=${DEFINES}' -DNP='${NP}' -DPLEV='${PLEV}' -I'${CONFIGH_PATH}

make $1
