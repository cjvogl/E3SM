#!/bin/bash

if [[ $# -lt 1 ]]; then
  echo "ERROR: specify a tstep_type number"
  exit -1
fi

if [[ $1 == 5 ]]; then
  TSTEP=(1.0 0.5 0.375 0.25 0.1875 0.125 0.0625 0.03125 0.015625 0.0078125 0.00390625 0.000390625)
  NMAX=(300 600 800 1200 1600 2400 4800 9600 19200 38400 76800 768000)
else
  TSTEP=(1.0 0.5 0.375 0.25 0.1875 0.125 0.0625 0.03125 0.015625 0.0078125 0.00390625)
  NMAX=(300 600 800 1200 1600 2400 4800 9600 19200 38400 76800)
fi

for n in ${!TSTEP[*]}; do
  if [[ $2 == "-noHV" ]]; then
    ./submit_gravitywave.py ./theta tsteptype $1 tstep ${TSTEP[$n]} nmax ${NMAX[$n]} nu 0.0
  else
    ./submit_gravitywave.py ./theta tsteptype $1 tstep ${TSTEP[$n]} nmax ${NMAX[$n]}
  fi
done
