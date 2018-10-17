#!/bin/bash

if [[ $# -lt 1 ]]; then
  echo "USAGE: ./run_convergence_test.sh <tstep_type value> [restart] [DCMIP X] [# levels] [rtol value] [-noHV]"
  exit -1
fi

TSTEP=(300 270 240 216 200 192 180 160 150 135 120 100 50 20 10)

for n in ${!TSTEP[*]}; do
  TYPEARG="tsteptype $1"
  if [[ $# -gt 1 ]]; then
    ARGS="ndays 1 restart $2"
  fi
  if [[ $# -gt 2 ]]; then
    VAR=$(bc <<< "scale=2;${TSTEP[$n]}/$3")
    TSTEPARG="$ARGS tstep $VAR"
    ARGS="$ARGS dcmip4_X $3"
  else
    TSTEPARG="tstep ${TSTEP[$n]}"
  fi
  EXE="theta-l-nlev"
  if [[ $# -gt 3 ]]; then
    EXE="${EXE}$4"
  else
    EXE="${EXE}30"
  fi
  if [[ $# -gt 4 ]]; then
    ARGS="$ARGS rtol $5"
  fi
  if [[ $4 == "-noHV" ]]; then
    ARGS="$ARGS nu 0.0"
  fi

  ./submit_baroclinicinstability.py $EXE $TYPEARG $TSTEPARG $ARGS
done
