#!/bin/bash

if [[ $# -lt 1 ]]; then
  echo "USAGE: ./run_convergence_test.sh <tstep_type value> [# levels] [rtol value] [-noHV]"
  exit -1
fi

# default levels is 30
if [[ $# -gt 2 ]]; then
  LEVELS=$2
else
  LEVELS=30
fi

TSTEP=(6.40 6.00 4.32 4.00 3.60 3.20 3.00 2.40 2.16 2.00 1.00)

for n in ${!TSTEP[*]}; do
  ARGS="tsteptype $1 tstep ${TSTEP[$n]}"
  if [[ $# -gt 2 ]]; then
    ARGS="$ARGS rtol $3"
  fi
  if [[ $4 == "-noHV" ]]; then
    ARGS="$ARGS nu 0.0"
  fi

  ./submit_baroclinicinstability.py theta-l-nlev$LEVELS $ARGS
done
