#!/bin/bash

WORKDIR=$S/homme/gravity_wave

if [[ $HOSTNAME == "cab"* ]]; then
  SYSTEM=cab
elif [[ $HOSTNAME == "quartz"* ]]; then
  SYSTEM=quartz
else
  echo "Script not implemented for this system"
  exit -1
fi

./cmake_homme_lc.sh 8 4 20 0 $1 $2

if [[ $? == 0 ]]; then
  mkdir -p $WORKDIR
  rm -f $WORKDIR/theta $WORKDIR/submit_gravitywave.py
  rm -f $WORKDIR/run_convergence_test.sh $WORKDIR/plot_convergence_test.py $WORKDIR/plot_performance.py
  cp build_$SYSTEM/src/theta/theta $WORKDIR/
  ln -s $PWD/dcmip_tests/dcmip2012_test3.1_nh_gravity_waves/theta/submit_gravitywave_quartz.py $WORKDIR/submit_gravitywave.py
  ln -s $PWD/dcmip_tests/dcmip2012_test3.1_nh_gravity_waves/theta/run_convergence_test_lc.sh $WORKDIR/run_convergence_test.sh
  ln -s $PWD/dcmip_tests/dcmip2012_test3.1_nh_gravity_waves/theta/plot_convergence_test_lc.py $WORKDIR/plot_convergence_test.py
  ln -s $PWD/dcmip_tests/dcmip2012_test3.1_nh_gravity_waves/theta/plot_performance_lc.py $WORKDIR/plot_performance.py
fi
