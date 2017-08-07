#!/bin/bash

WORKDIR=$S/homme/jw_baroclinic

if [[ $HOSTNAME == "cab"* ]]; then
  SYSTEM=cab
elif [[ $HOSTNAME == "quartz"* ]]; then
  SYSTEM=quartz
else
  echo "Script not implemented for this system"
  exit -1
fi

./cmake_homme_lc.sh 8 4 26 4

if [[ $? == 0 ]]; then
  mkdir -p $WORKDIR
  rm -f $WORKDIR/theta $WORKDIR/baro.job $WORKDIR/jw_baroclini.nl $WORKDIR/*.ascii
  cp build_$SYSTEM/src/theta/theta $WORKDIR/
  cp test/jw_baroclinic/baro.job $WORKDIR/
  cp test/jw_baroclinic/jw_baroclinic.nl $WORKDIR/
  cp test/vcoord/*.ascii $WORKDIR/
fi
