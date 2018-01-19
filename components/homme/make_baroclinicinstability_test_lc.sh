#!/bin/bash

WORKDIR=$S/homme/baroclinicinstability_test
NLEV=30

if [[ $HOSTNAME == "cab"* ]]; then
  SYSTEM=cab
elif [[ $HOSTNAME == "quartz"* ]]; then
  SYSTEM=quartz
else
  echo "Script not implemented for this system"
  exit -1
fi

./cmake_homme_lc.sh 8 4 $NLEV 0 $1

if [[ $? == 0 ]]; then
  mkdir -p $WORKDIR
  rm -f $WORKDIR/theta-nlev$NLEV $WORKDIR/submit_baroclinicinstability.py $WORKDIR/find_largest_timestep.py
  cp build_$SYSTEM/src/theta/theta $WORKDIR/theta-nlev$NLEV
  cp build_$SYSTEM/dcmip_tests/dcmip2012_test4.1_baroclinic_instability/vcoord/camm-30.ascii $WORKDIR/
  cp build_$SYSTEM/dcmip_tests/dcmip2012_test4.1_baroclinic_instability/vcoord/cami-30.ascii $WORKDIR/
  ln -s $PWD/dcmip_tests/dcmip2012_test4.1_baroclinic_instability/theta/submit_baroclinicinstability_lc.py $WORKDIR/submit_baroclinicinstability.py
  ln -s $PWD/dcmip_tests/dcmip2012_test4.1_baroclinic_instability/theta/find_largest_timestep_lc.py $WORKDIR/find_largest_timestep.py
fi
