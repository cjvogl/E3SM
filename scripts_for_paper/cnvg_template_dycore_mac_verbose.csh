#====================================================================
# Create and set up new case. No need to build the model. 
#====================================================================
source ${case_setup_script}

cd $CASEROOT

echo ./xmlchange -file env_build.xml -id BUILD_COMPLETE  -val 'TRUE'
./xmlchange -file env_build.xml -id BUILD_COMPLETE  -val 'TRUE'

#-----------------------------------------
# Runtime options: edit env_run.xml
#-----------------------------------------
cd $CASEROOT

@ nresub = $ncycle - 1

echo ./xmlchange  -file env_run.xml -id  STOP_N       -val $nlen
./xmlchange  -file env_run.xml -id  STOP_N       -val $nlen
echo ./xmlchange  -file env_run.xml -id  STOP_OPTION  -val $stop_o 
./xmlchange  -file env_run.xml -id  STOP_OPTION  -val $stop_o 
echo ./xmlchange  -file env_run.xml -id  REST_N       -val $nlen
./xmlchange  -file env_run.xml -id  REST_N       -val $nlen
echo ./xmlchange  -file env_run.xml -id  REST_OPTION  -val 'ndays' 
./xmlchange  -file env_run.xml -id  REST_OPTION  -val 'ndays' 
echo ./xmlchange  -file env_run.xml -id  RESUBMIT     -val $nresub 
./xmlchange  -file env_run.xml -id  RESUBMIT     -val $nresub 
echo ./xmlchange  -file env_run.xml -id  DOUT_S       -val 'FALSE'
./xmlchange  -file env_run.xml -id  DOUT_S       -val 'FALSE'
echo ./xmlchange  -file env_run.xml -id  DOUT_L_MS    -val 'FALSE'
./xmlchange  -file env_run.xml -id  DOUT_L_MS    -val 'FALSE'

set ncpl = `echo 86400 / $dtime | bc`

echo ./xmlchange  -file env_run.xml -id  ATM_NCPL          -val $ncpl
./xmlchange  -file env_run.xml -id  ATM_NCPL          -val $ncpl
echo ./xmlchange  -file env_run.xml -id  CAM_NAMELIST_OPTS -val dtime=$dtime 
./xmlchange  -file env_run.xml -id  CAM_NAMELIST_OPTS -val dtime=$dtime 
echo ./xmlchange  -file env_run.xml -id  CLM_NAMELIST_OPTS -val dtime=$dtime 
./xmlchange  -file env_run.xml -id  CLM_NAMELIST_OPTS -val dtime=$dtime 

#--------------------
# Namelist variables
#--------------------
cat <<EOF >> user_nl_cam
 avgflag_pertape    = 'I'
 nhtfrq             = $nhtfrq,
 mfilt              = $mfilt,
 ndens              = 1,
 empty_htapes       = .true.,
 fincl1             = 'PS','U','V','T','Q','CLDLIQ','RKZ_fac','RKZ_f','RKZ_RH','RKZ_Av','RKZ_Al','RKZ_AT','RKZ_qlneg_lm4','RKZ_qvneg_lm4','RKZ_qme','RKZ_term_A','RKZ_term_B','RKZ_term_C'
 ncdata             = '$atm_init'
 deep_scheme        = 'off',
 shallow_scheme     = 'off',
 l_tracer_aero      = .false.
 l_vdiff            = .false.
 l_rayleigh         = .false.
 l_gw_drag          = .false.
 l_ac_energy_chk    = .true.
 l_bc_energy_fix    = .true.
 l_dry_adj          = .false.
 l_st_mac           = .true.
 l_st_mic           = .false.
 l_rad              = .false.
EOF
cat $nl_file >> user_nl_cam

#============
# Run model 
#============
cd $CASEROOT

if ( $postCIME <= 2 ) then
   ./$CASE.submit
else
   ./case.submit
endif

date
