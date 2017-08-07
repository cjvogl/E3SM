#!/bin/bash

WORKDIR=$S/homme/dcmip_test1

if [[ $HOSTNAME == "quartz"* ]]; then
  SYSTEM=quartz
else
  echo "Script not implemented for this system"
  exit -1
fi

./cmake_homme_lc.sh 4 4 30 0

if [[ $? == 0 ]]; then
  mkdir -p $WORKDIR
  cp build_$SYSTEM/src/theta/theta $WORKDIR/theta_$SYSTEM
  cp dcmip_tests/dcmip2016_test1_baroclinic_wave/theta/namelist-nh.nl $WORKDIR
  cp dcmip_tests/dcmip2016_test1_baroclinic_wave/theta/jobscript-lc.sh $WORKDIR
fi
 
