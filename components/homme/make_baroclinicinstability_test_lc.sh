#!/bin/bash

WORKDIR=$S/homme/baroclinicinstability_test_noheating
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
  rm -f $WORKDIR/theta-l-nlev$NLEV $WORKDIR/submit_baroclinicinstability.py
  rm -f $WORKDIR/run_convergence_test.sh $WORKDIR/plot_convergence_test.py
  rm -f $WORKDIR/plot_efficiency.py $WORKDIR/plot_energy_test.py $WORKDIR/plot_energy_efficiency.py
  rm -f $WORKDIR/dictionaries.py
  cp build_$SYSTEM/src/theta-l/theta-l $WORKDIR/theta-l-nlev$NLEV
  cp build_$SYSTEM/dcmip_tests/dcmip2012_test4.1_baroclinic_instability/vcoord/camm-30.ascii $WORKDIR/
  cp build_$SYSTEM/dcmip_tests/dcmip2012_test4.1_baroclinic_instability/vcoord/cami-30.ascii $WORKDIR/
  ln -s $PWD/dcmip_tests/dcmip2012_test4.1_baroclinic_instability/theta-l/submit_baroclinicinstability_lc.py $WORKDIR/submit_baroclinicinstability.py
  ln -s $PWD/dcmip_tests/dcmip2012_test4.1_baroclinic_instability/theta-l/run_convergence_test_lc.sh $WORKDIR/run_convergence_test.sh
  ln -s $PWD/plot_scripts/plot_convergence_test.py $WORKDIR/plot_convergence_test.py
  ln -s $PWD/plot_scripts/plot_efficiency.py $WORKDIR/plot_efficiency.py
  ln -s $PWD/plot_scripts/plot_energy_test.py $WORKDIR/plot_energy_test.py
  ln -s $PWD/plot_scripts/plot_energy_efficiency.py $WORKDIR/plot_energy_efficiency.py
  ln -s $PWD/plot_scripts/dictionaries.py $WORKDIR/dictionaries.py
fi
