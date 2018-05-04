#!/bin/bash

if [[ $# -lt 1 ]]; then
  echo "USAGE: ./run_convergence_test.sh <tstep_type value> [rtol value] [imex splitting] [-noHV]"
  exit -1
fi

if [[ $1 == 5 ]]; then
  TSTEP=(600 400 200 100 10)
else
  TSTEP=(640 600 432 400 360 320 300 240 216 200 100)
fi

for n in ${!TSTEP[*]}; do
  ARGS="tsteptype $1 tstep ${TSTEP[$n]}"
  if [[ $# -gt 1 ]]; then
    ARGS="$ARGS rtol $2"
  fi
  if [[ $# -gt 2 ]]; then
    ARGS="$ARGS splitting $3"
  fi
  if [[ $4 == "-noHV" ]]; then
    ARGS="$ARGS nu 0.0"
  fi
  if [[ $1 == 5 ]]; then
    ARGS="$ARGS hydrostatic true"
  fi

  ./submit_baroclinicinstability.py theta-l-nlev30 $ARGS
done
