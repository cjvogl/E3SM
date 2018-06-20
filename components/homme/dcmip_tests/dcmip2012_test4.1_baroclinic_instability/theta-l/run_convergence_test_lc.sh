#!/bin/bash

if [[ $# -lt 1 ]]; then
  echo "USAGE: ./run_convergence_test.sh <tstep_type value> [DCMIP X] [# levels] [rtol value] [-noHV]"
  exit -1
fi

TSTEP=(300 270 240 216 200 192 180 160 150 135 120 100 50 20 10)

for n in ${!TSTEP[*]}; do
  ARGS="tsteptype $1"
  TSTEPCURRENT=${TSTEP[$n]}
  if [[ $# -gt 1 ]]; then
    VAR=$(bc <<< "scale=2;${TSTEP[$n]}/$2")
    ARGS="$ARGS tstep $VAR dcmip4_X $2"
  else
    ARGS="$ARGS tstep ${TSTEP[$n]} dcmip4_X 1"
  fi
  EXE="theta-l-nlev"
  if [[ $# -gt 2 ]]; then
    EXE="${EXE}$3"
  else
    EXE="${EXE}30"
  fi
  if [[ $# -gt 3 ]]; then
    ARGS="$ARGS rtol $4"
  fi
  if [[ $4 == "-noHV" ]]; then
    ARGS="$ARGS nu 0.0"
  fi

  ./submit_baroclinicinstability.py $EXE $ARGS
done
