#====================================================================
# create new case
#====================================================================
rm -rf $CASEROOT

if ( $postCIME == 0 ) then
    cd  $CCSMROOT/scripts
else
    cd  $CCSMROOT/cime/scripts
endif

./create_newcase -case $CASEROOT -mach $MACH \
                 -res $RESOLUTION -compset $COMPSET -compiler $COMPILER #-v

#====================================================================
# set up case
#====================================================================
cd $CASEROOT

echo ./xmlchange -file env_run.xml   -id RUNDIR  -val $RUNDIR
./xmlchange -file env_run.xml   -id RUNDIR  -val $RUNDIR
echo ./xmlchange -file env_build.xml -id EXEROOT -val $EXEDIR
./xmlchange -file env_build.xml -id EXEROOT -val $EXEDIR

echo ./xmlchange -file env_run.xml   -id DIN_LOC_ROOT          -val $CSMDATA
./xmlchange -file env_run.xml   -id DIN_LOC_ROOT          -val $CSMDATA
echo ./xmlchange -file env_run.xml   -id DIN_LOC_ROOT_CLMFORC  -val $CSMDATA
./xmlchange -file env_run.xml   -id DIN_LOC_ROOT_CLMFORC  -val $CSMDATA

if ( $QUEUE == "regular" ) then
     echo ./xmlchange -file env_run.xml -id USER_REQUESTED_QUEUE -val "regular" 
     ./xmlchange -file env_run.xml -id USER_REQUESTED_QUEUE -val "regular" 
     echo ./xmlchange -file env_run.xml -id USER_REQUESTED_WALLTIME -val "10:00:00"
     ./xmlchange -file env_run.xml -id USER_REQUESTED_WALLTIME -val "10:00:00"
else if ( $QUEUE == "debug" ) then
     if ( $MACH == "quartz" ) then
          echo ./xmlchange -file env_run.xml -id USER_REQUESTED_QUEUE -val "pdebug"
          ./xmlchange -file env_run.xml -id USER_REQUESTED_QUEUE -val "pdebug"
     else
          echo ./xmlchange -file env_run.xml -id USER_REQUESTED_QUEUE -val "debug"
          ./xmlchange -file env_run.xml -id USER_REQUESTED_QUEUE -val "debug"
     endif
     echo ./xmlchange -file env_run.xml -id USER_REQUESTED_WALLTIME -val "00:30:00"
     ./xmlchange -file env_run.xml -id USER_REQUESTED_WALLTIME -val "00:30:00"
endif

# set the PIO_TYPE (needed for running on quartz at LLNL)
if ($MACH == "quartz") then
    echo ./xmlchange -file env_run.xml -id PIO_TYPENAME -val 'netcdf'
    ./xmlchange -file env_run.xml -id PIO_TYPENAME -val 'netcdf'
endif

# Specify PE layout and # of instances

@ nproc = $NTASKS_PER_INST * $NINST

echo ./xmlchange -file env_mach_pes.xml -id NTASKS_ATM -val $nproc
./xmlchange -file env_mach_pes.xml -id NTASKS_ATM -val $nproc
echo ./xmlchange -file env_mach_pes.xml -id NTHRDS_ATM -val $NTHRDS
./xmlchange -file env_mach_pes.xml -id NTHRDS_ATM -val $NTHRDS
echo ./xmlchange -file env_mach_pes.xml -id ROOTPE_ATM -val '0'
./xmlchange -file env_mach_pes.xml -id ROOTPE_ATM -val '0'
echo ./xmlchange -file env_mach_pes.xml -id NINST_ATM -val $NINST
./xmlchange -file env_mach_pes.xml -id NINST_ATM -val $NINST
echo ./xmlchange -file env_mach_pes.xml -id NINST_ATM_LAYOUT -val 'concurrent'
./xmlchange -file env_mach_pes.xml -id NINST_ATM_LAYOUT -val 'concurrent'

echo ./xmlchange -file env_mach_pes.xml -id NTASKS_LND -val $nproc
./xmlchange -file env_mach_pes.xml -id NTASKS_LND -val $nproc
echo ./xmlchange -file env_mach_pes.xml -id NTHRDS_LND -val $NTHRDS
./xmlchange -file env_mach_pes.xml -id NTHRDS_LND -val $NTHRDS
echo ./xmlchange -file env_mach_pes.xml -id ROOTPE_LND -val '0'
./xmlchange -file env_mach_pes.xml -id ROOTPE_LND -val '0'
echo ./xmlchange -file env_mach_pes.xml -id NINST_LND  -val $NINST
./xmlchange -file env_mach_pes.xml -id NINST_LND  -val $NINST
echo ./xmlchange -file env_mach_pes.xml -id NINST_LND_LAYOUT -val 'concurrent'
./xmlchange -file env_mach_pes.xml -id NINST_LND_LAYOUT -val 'concurrent'

echo ./xmlchange -file env_mach_pes.xml -id NTASKS_ROF -val $nproc
./xmlchange -file env_mach_pes.xml -id NTASKS_ROF -val $nproc
echo ./xmlchange -file env_mach_pes.xml -id NTHRDS_ROF -val $NTHRDS
./xmlchange -file env_mach_pes.xml -id NTHRDS_ROF -val $NTHRDS
echo ./xmlchange -file env_mach_pes.xml -id ROOTPE_ROF -val '0'
./xmlchange -file env_mach_pes.xml -id ROOTPE_ROF -val '0'
echo ./xmlchange -file env_mach_pes.xml -id NINST_ROF  -val $NINST
./xmlchange -file env_mach_pes.xml -id NINST_ROF  -val $NINST
echo ./xmlchange -file env_mach_pes.xml -id NINST_ROF_LAYOUT -val 'concurrent'
./xmlchange -file env_mach_pes.xml -id NINST_ROF_LAYOUT -val 'concurrent'

echo ./xmlchange -file env_mach_pes.xml -id NTASKS_ICE -val $nproc
./xmlchange -file env_mach_pes.xml -id NTASKS_ICE -val $nproc
echo ./xmlchange -file env_mach_pes.xml -id NTHRDS_ICE -val $NTHRDS 
./xmlchange -file env_mach_pes.xml -id NTHRDS_ICE -val $NTHRDS 
echo ./xmlchange -file env_mach_pes.xml -id ROOTPE_ICE -val '0'
./xmlchange -file env_mach_pes.xml -id ROOTPE_ICE -val '0'
echo ./xmlchange -file env_mach_pes.xml -id NINST_ICE -val $NINST
./xmlchange -file env_mach_pes.xml -id NINST_ICE -val $NINST
echo ./xmlchange -file env_mach_pes.xml -id NINST_ICE_LAYOUT -val 'concurrent'
./xmlchange -file env_mach_pes.xml -id NINST_ICE_LAYOUT -val 'concurrent'

echo ./xmlchange -file env_mach_pes.xml -id NTASKS_OCN -val $nproc
./xmlchange -file env_mach_pes.xml -id NTASKS_OCN -val $nproc
echo ./xmlchange -file env_mach_pes.xml -id NTHRDS_OCN -val $NTHRDS 
./xmlchange -file env_mach_pes.xml -id NTHRDS_OCN -val $NTHRDS 
echo ./xmlchange -file env_mach_pes.xml -id ROOTPE_OCN -val '0'
./xmlchange -file env_mach_pes.xml -id ROOTPE_OCN -val '0'
echo ./xmlchange -file env_mach_pes.xml -id NINST_OCN -val $NINST
./xmlchange -file env_mach_pes.xml -id NINST_OCN -val $NINST
echo ./xmlchange -file env_mach_pes.xml -id NINST_OCN_LAYOUT -val 'concurrent'
./xmlchange -file env_mach_pes.xml -id NINST_OCN_LAYOUT -val 'concurrent'

# GLC and WAV are stub components
echo ./xmlchange -file env_mach_pes.xml -id NTASKS_GLC -val $nproc
./xmlchange -file env_mach_pes.xml -id NTASKS_GLC -val $nproc
echo ./xmlchange -file env_mach_pes.xml -id NTHRDS_GLC -val $NTHRDS 
./xmlchange -file env_mach_pes.xml -id NTHRDS_GLC -val $NTHRDS 
echo ./xmlchange -file env_mach_pes.xml -id ROOTPE_GLC -val '0'
./xmlchange -file env_mach_pes.xml -id ROOTPE_GLC -val '0'
echo ./xmlchange -file env_mach_pes.xml -id NINST_GLC -val '1'
./xmlchange -file env_mach_pes.xml -id NINST_GLC -val '1'
echo ./xmlchange -file env_mach_pes.xml -id NINST_GLC_LAYOUT -val 'concurrent'
./xmlchange -file env_mach_pes.xml -id NINST_GLC_LAYOUT -val 'concurrent'

echo ./xmlchange -file env_mach_pes.xml -id NTASKS_WAV -val $nproc
./xmlchange -file env_mach_pes.xml -id NTASKS_WAV -val $nproc
echo ./xmlchange -file env_mach_pes.xml -id NTHRDS_WAV -val $NTHRDS 
./xmlchange -file env_mach_pes.xml -id NTHRDS_WAV -val $NTHRDS 
echo ./xmlchange -file env_mach_pes.xml -id ROOTPE_WAV -val '0'
./xmlchange -file env_mach_pes.xml -id ROOTPE_WAV -val '0'
echo ./xmlchange -file env_mach_pes.xml -id NINST_WAV -val '1'
./xmlchange -file env_mach_pes.xml -id NINST_WAV -val '1'
echo ./xmlchange -file env_mach_pes.xml -id NINST_WAV_LAYOUT -val 'concurrent'
./xmlchange -file env_mach_pes.xml -id NINST_WAV_LAYOUT -val 'concurrent'

echo ./xmlchange -file env_mach_pes.xml -id NTASKS_CPL -val $nproc
./xmlchange -file env_mach_pes.xml -id NTASKS_CPL -val $nproc
echo ./xmlchange -file env_mach_pes.xml -id NTHRDS_CPL -val $NTHRDS 
./xmlchange -file env_mach_pes.xml -id NTHRDS_CPL -val $NTHRDS 
echo ./xmlchange -file env_mach_pes.xml -id ROOTPE_CPL -val '0'
./xmlchange -file env_mach_pes.xml -id ROOTPE_CPL -val '0'

if ( $postCIME <= 2 ) then
   ./cesm_setup
else
   ./case.setup
endif
