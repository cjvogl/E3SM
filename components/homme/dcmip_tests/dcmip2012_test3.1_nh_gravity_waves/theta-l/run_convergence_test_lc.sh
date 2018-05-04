#!/bin/bash

if [[ $# -lt 1 ]]; then
  echo "USAGE: ./run_convergence_test.sh <tstep_type value> [rtol value] [imex splitting] [-noHV]"
  exit -1
fi

if [[ $1 == 5 ]]; then
  TSTEP=(4.0 3.0 2.0 1.0 0.5 0.375 0.25 0.1875 0.125 0.0625 0.03125 0.015625 0.0078125 0.00390625 0.000390625)
  NMAX=(75 100 150 300 600 800 1200 1600 2400 4800 9600 19200 38400 76800 768000)
else
  TSTEP=(4.0 3.0 2.0 1.0 0.5 0.375 0.25 0.1875 0.125 0.0625 0.03125 0.015625 0.0078125 0.00390625)
  NMAX=(75 100 150 300 600 800 1200 1600 2400 4800 9600 19200 38400 76800)
fi

for n in ${!TSTEP[*]}; do
  ARGS="tsteptype $1 tstep ${TSTEP[$n]} nmax ${NMAX[$n]}"
  if [[ $# -gt 1 ]]; then
    ARGS="$ARGS rtol $2"
  fi
  if [[ $# -gt 2 ]]; then
    ARGS="$ARGS splitting $3"
  fi
  if [[ $4 == "-noHV" ]]; then
    ARGS="$ARGS nu 0.0"
  fi

  ./submit_gravitywave.py theta-l-nlev20 $ARGS
done
