#!/bin/bash

WORKDIR=$S/homme/baroclinicwave_test
NLEV=30

if [[ $HOSTNAME == "cab"* ]]; then
  SYSTEM=cab
elif [[ $HOSTNAME == "quartz"* ]]; then
  SYSTEM=quartz
else
  echo "Script not implemented for this system"
  exit -1
fi

./cmake_homme_lc.sh 8 4 $NLEV 5 $1

if [[ $? == 0 ]]; then
  mkdir -p $WORKDIR
  rm -f $WORKDIR/theta-nlev$NLEV $WORKDIR/submit_baroclinicwave.py
#  rm -f $WORKDIR/run_convergence_test.sh $WORKDIR/plot_convergence_test.py $WORKDIR/plot_performance.py
  cp build_$SYSTEM/src/theta/theta $WORKDIR/theta-nlev$NLEV
  cp build_$SYSTEM/dcmip_tests/dcmip2016_test1_baroclinic_wave/vcoord/camm-30.ascii $WORKDIR/
  cp build_$SYSTEM/dcmip_tests/dcmip2016_test1_baroclinic_wave/vcoord/cami-30.ascii $WORKDIR/
  ln -s $PWD/dcmip_tests/dcmip2016_test1_baroclinic_wave/theta/submit_baroclinicwave_lc.py $WORKDIR/submit_baroclinicwave.py
#  ln -s $PWD/dcmip_tests/dcmip2012_test3.1_nh_gravity_waves/theta/run_convergence_test_lc.sh $WORKDIR/run_convergence_test.sh
#  ln -s $PWD/dcmip_tests/dcmip2012_test3.1_nh_gravity_waves/theta/plot_convergence_test_lc.py $WORKDIR/plot_convergence_test.py
#  ln -s $PWD/dcmip_tests/dcmip2012_test3.1_nh_gravity_waves/theta/plot_performance_lc.py $WORKDIR/plot_performance.py
fi
