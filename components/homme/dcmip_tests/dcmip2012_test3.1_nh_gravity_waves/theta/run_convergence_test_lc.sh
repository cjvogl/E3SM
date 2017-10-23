#!/bin/bash

if [[ $HOSTNAME == "cab"* ]]; then
  SYSTEM=cab
elif [[ $HOSTNAME == "quartz"* ]]; then
  SYSTEM=quartz
else
  echo "ERROR: Script not implemented for this system"
  exit -1
fi

if [[ $# -lt 1 ]]; then
  echo "ERROR: specify a tstep_type number"
  exit -1
fi

if [[ $1 == 5 ]]; then
  TSTEP=(1.0 0.5 0.375 0.25 0.1875 0.125 0.0625 0.03125 0.015625 0.0015625)
  NMAX=(600 1200 1600 2400 3200 4800 9600 19200 38400 384000)
else
  TSTEP=(1.0 0.5 0.375 0.25 0.1875 0.125 0.0625 0.03125 0.015625)
  NMAX=(600 1200 1600 2400 3200 4800 9600 19200 38400)
fi

for n in ${!TSTEP[*]}; do
  if [[ $2 == "-noHV" ]]; then
    ./submit_gravitywave.py ./theta tsteptype $1 tstep ${TSTEP[$n]} nmax ${NMAX[$n]} nu 0.0
  else
    ./submit_gravitywave.py ./theta tsteptype $1 tstep ${TSTEP[$n]} nmax ${NMAX[$n]}
  fi
done
