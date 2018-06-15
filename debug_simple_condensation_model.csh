#!/bin/csh
date

##############################################################################################
# This is the driver scripts for the time step convergence test simulations which
#  1. gets the E3SM code from GitHub;
#  2. compiles the model (just once);
#  3. generates one script for each simulation, then runs those scripts to create cases 
#     submit jobs. 
# You can turn on or off these steps using the variables below. 1 = on; 0 = off.
# Note that for step 1 you have two options: clone the code into a new local repo using
# "git clone" (switch clone_code) or update an existing local repo using "git pull"
# (switch update_code). 
#
# Script history:
# Original and revised versions: Hui Wan (PNNL), 2014-2017.
# Clean-up for the SciDAC Convergence project repo: Hui Wan (PNNL), 2017-12-01 
#
# Contact: Hui.Wan@pnnl.gov
##############################################################################################
set clone_code    = 0
set update_code   = 0

set compile_model = 0
set run_model     = 1

#-----------------------------------------------------------------------------------------
# The variable $taskname is used to organize the exe/run/case directories (and is NOT used 
# for other purposes). Because the convergence tests typically involve a large number of 
# simulations that can share the same executable, I (Hui) decided to put the exe/run/case
# directories all in one place in a large patition of the system's disk space that 
# does not automatically purge files after a certain period of time. I create cases 
# with names starting with "compile_" to compile the model, then create separate cases
# with different names (that never start with "compile_" for conducting simulations to 
# test solution convergence. Also, unlike the typical workflow described in the CESM tutorials,
# usually do not let the model generate  restart files for these simulations because the 
# runs  are very short) and I leave the output files in the run directory (i.e., short-term 
# archiving is turned off).
#
# The root directory for this "task" which can contain many sets of convergence tests
# using different model configurations is 
#   $PTMP/$taskname/
# where $PTMP is defined further below. Under $PTMP/$taskname/, 
# three subdirectories are created by this script:
#
#  - exe/
#  - cases/
#  - run/
#
# Under exe/, there might be only one subdirectory called compile_... if the E3SM model 
# is compiled # for only one compset/resolution/PE-layout; there might also be multiple
# subdirectories if the model is compiled for different configurations. 
# These subdirectories contain the corresponding model executables.
#
# Under cases/ and run/, you will find sub-directories for both the compilation and 
# the test simulations, althought under run/compile_* there are no model output
# because the compile_* cases are used for model compilation only.
#
# If you want to put the exe/ directory somewhere else, please change the value of EXELOC.
# 
set taskname = cnvg_condensation

# The variable $testconfig is used to pick the appropriate template file that contains certain
# run-time options (amongh those are namelist variables that specify the contents and 
# frequency of model output, and switches that turn on or off certain parameterizations)
#
#set testconfig   = 'dycore_only'
set testconfig   = 'dycore_mac'    # dycore + mac

# Compset, machine, resolution, etc.
#
setenv COMPSET    FC5AQUAP

setenv MACH       edison
setenv CESM_EMAIL vogl2@llnl.gov
setenv RESOLUTION ne30_ne30

####################################################################
# Clone or update code
####################################################################
setenv CCSMTAG e3sm 
setenv CCSMROOT $HOME/workspace/${CCSMTAG}
setenv BRANCH 'huiwanpnnl/atm/condensation'

# CIME was upgraded from version 2 to version 5 in early 2017 
# which led to name and/or path changes for a number of CIME 
# scripts. The variable postCIME was added here for backward compatibility.
# For earlier versions of the model, please set postCIME to 2;
# For codes from mid/late 2017, please set postCIME to 5.
#
setenv postCIME 5

# Clone code repo

if ($clone_code > 0) then

   mkdir -p $CCSMROOT:h
   cd $CCSMROOT:h
   git clone git@github.com:ACME-Climate/ACME.git $CCSMTAG

   cd $CCSMTAG
   git submodule update --init
   git checkout $BRANCH

endif

# Or - update the existing repo

if ($update_code > 0) then

   cd $CCSMROOT
   git pull

endif

#############################################################################
# Conduct simulations with different model configurations
#----------------------------------------------------------------------------
# For test simulations, the casenames have the following pattern:
#
#  ${group}_${RESOLUTION}_${COMPSET}_${COMPILER}_${MACH}_"DT"`printf "%04d" ${dtime}`"_"${irstring}_${testconfig}
#
# where ${group} is one of the elements in $groupList defined immediately below,
# and ${testconfig} was defined earlier in this script.
#
# Note that you can put additional namelist settings into the file
#
#   ./namelist_files/cam_nl_$group
#
# to further distinguish different simulations. The differences between 
# $testconfig and $group are: 
#
#   -  for every $testconfig you need a file cnvg_template_${testconfig}.csh
#      which will become part of the shell script that runs a test simulationc;
#      the run-time settings in there can contain both xml changes and 
#      namelist variable specifications in user_nl_cam etc.
#
#   -  for every $group you need a file ./namelist_files/cam_nl_$group
#      which can contain only namelist variable specifications. Those will 
#      end up in user_nl_cam.
#
# For the namelist variables, I personally use $testconfig to turn on or off
# certain parameterizations and adjust the output field list if needed. 
# For example, I use different $testconfig for test simulations with dycore only
# or with dycore + macrophysics; When I want to further distinguish different
# configurations of, say, the simple condensation model (e.g., use different
# options for a term on the RHS of the condensation equation), those disctions
# are made using $group.
#
# If you do not have anything extra to specify for $group and do not create
# the file ./namelist_files/cam_nl_$group, the script will complain about that,
# but no harm will be done.
#
# DEBUG
#set groupList = ("RKZ_A1_B0_C0_lmt4")

#set groupList = ("RKZ_A1_B1_C0_lmt4")

#set groupList = ("RKZ_A1_B0_C1_ql7_lmt4")
#set groupList = ("RKZ_A1_B0_C2_ql17_lmt4")

#set groupList = ("RKZ_A1_B1_C1_ql7_lmt4")

#set groupList = ("RKZ_A1_B1_C2_ql17_lmt4")
#set groupList = ("RKZ_A1_B1_C2_ql19_lmt4")

set groupList = ("RKZ_P1_lmt4")
#set groupList = ("RKZ_P2_lmt4")
#set groupList = ("RKZ_P1_lmt4_adjIC")
#set groupList = ("RKZ_P2_lmt4_adjIC")

# If multiple groups are specified above, you can use igS (group index start)
# and igE (group index end) to run a subset in one particulation execution of 
# this script.
#
set ngroups = $#groupList
set igS = 1
set igE = $ngroups

#-------------------------------------------------------------------------------
# Start and end indicies of ensemble members. Ensemble members differ in the 
# initial conditions. For each ensemble member, the same initial conditions 
# are used to conduct multiple simulations with different time step sizes
# specified by $dtimeList further below.
#-------------------------------------------------------------------------------
set irS = 1   # ensemble member: start index
# DEBUG
set irE = 1   # ensemble member: end   index
#set irE = 6   # ensemble member: end   index

#-------------------------------------------------------------------------------
# List of time step sizes with with simulations will be conducted (unit: second)
#-------------------------------------------------------------------------------
# DEBUG
set dtimeList = (1800)
set ndtime = $#dtimeList

#-------------------------------------------------------------------------------
# Compile and run model in debug mode?
#-------------------------------------------------------------------------------
# DEBUG
#set debug = 'TRUE'
set debug = 'FALSE'

#====================================================================
# Paths to source code and model input/output
#====================================================================
setenv NINST  1

if ($MACH == "titan") then

   setenv CESM_PROJ cli112
   setenv PROJECT $CESM_PROJ
   setenv COMPILER   pgi
   setenv COMPILER   intel
   setenv NTHRDS 1
   setenv NTASKS_PER_INST 512

   setenv CSMDATA  /lustre/atlas1/cli900/world-shared/cesm/inputdata
   setenv PTMP     /lustre/atlas/proj-shared/cli112/huiwan/$taskname
   setenv EXELOC   $PTMP/exe/
   set initDir = /lustre/atlas1/cli112/proj-shared/huiwan/cesm_input/FC5AV1C-02_init/

else if ($MACH == "eos") then

   setenv CESM_PROJ cli112
   setenv PROJECT $CESM_PROJ
   setenv COMPILER   pgi
   setenv COMPILER   intel
   setenv NTHRDS 4
   setenv NTASKS_PER_INST 64

   setenv CSMDATA  /lustre/atlas1/cli900/world-shared/cesm/inputdata
   setenv PTMP     /lustre/atlas/proj-shared/cli112/huiwan/$taskname
   setenv EXELOC   $PTMP/exe/
   set initDir = /lustre/atlas1/cli112/proj-shared/huiwan/cesm_input/FC5AV1C-02_init/ 

else if ($MACH == "corip1") then  # this is probably outdated.

   setenv CESM_PROJ  m1199
   setenv PROJECT $CESM_PROJ
   setenv COMPILER   intel
   setenv NTASKS_PER_INST 384

   setenv CSMDATA  /project/projectdirs/acme/inputdata/ 
   setenv PTMP   /global/cscratch1/sd/huiwan/huiwan/$taskname
   setenv EXELOC /global/cscratch1/sd/huiwan/huiwan/$taskname/exe
   set initDir = /global/cscratch1/sd/huiwan/cesm_input/acme_capt_init_FC5_weekly_sst/

else if ($MACH == "constance") then

   setenv CESM_PROJ uq_climate
   setenv PROJECT $CESM_PROJ
   setenv COMPILER   intel
   setenv NTHRDS 1
   setenv NTASKS_PER_INST 288

   setenv CSMDATA  /pic/projects/climate/csmdata/ 
   setenv PTMP     /pic/projects/uq_climate/wanh895/$taskname
   setenv EXELOC   $PTMP/exe/
   set initDir = /pic/projects/uq_climate/wanh895/acme_input/ne30_FC5_init/ 

else if ($MACH == "quartz") then

   setenv CESM_PROJ climalg
   setenv PROJECT $CESM_PROJ
   setenv COMPILER intel
   setenv NTHRDS 1
   setenv NTASKS_PER_INST 288

   setenv CSMDATA  /usr/gdata/climdat/ccsm3data/inputdata/
   setenv PTMP     /p/lscratchh/$USER/ACME/$taskname
   setenv EXELOC   $PTMP/exe
   set initDir = /p/lscratchh/$USER/acme_input/ne30_FC5_init/

else if ($MACH == "edison") then
   
   setenv CESM_PROJ m3089
   setenv PROJECT $CESM_PROJ
   setenv COMPILER intel
   setenv NTHRDS 1
   setenv NTASKS_PER_INST 288
   setenv QUEUE "debug"

   setenv CSMDATA /project/projectdirs/acme/inputdata/
   setenv PTMP $SCRATCH/e3sm/$taskname
   setenv EXELOC $PTMP/exe
   set initDir = /global/cscratch1/sd/zhan391/acme_input/ne30_FC5_init/

else

   echo Specify paths for MACH $MACH ! Abort.
   exit

endif
   mkdir -p $PTMP 

#----------------------------------------------------------------------------------
# We will compile the model just once, and use the same executable for all ensemble
# members. Personally, I prefer to create a separate "case" just for the compilation.

if ($NINST > 1) then
   set execase = compile_${CCSMTAG}_${COMPSET}_${RESOLUTION}_${MACH}_${COMPILER}_${NINST}x${NTASKS_PER_INST}bundle
else
   set execase = compile_${CCSMTAG}_${COMPSET}_${RESOLUTION}_${MACH}_${COMPILER}_${NTASKS_PER_INST}proc
endif

if ($debug == 'TRUE') then
   set execase = ${execase}_debug
endif

setenv EXEDIR ${EXELOC}/$execase

#----------------------------------------------------------------------------
# Later in this script, one "case" is created for every ensemble member (realization)
# of the reference simulations, trusted simulations, and test simulations. 
# We need to make sure that the same compile-time options (e.g., compset, 
# compiler, PE layout) are used for those "simulation cases" as well as 
# for the "compilation case". $case_setup_script contains the "create_case"
# command, specification of RUNDIR, EXEROOT, and PE layout, as well as 
# the "cesm_setup" command. This script is sourced (i.e., used like a
# Fortran subroutine) when we create new case for compilation and simulations.

set exp_setup_dir = `pwd`
set case_setup_script = ${exp_setup_dir}/create_and_setup_bundled_case.csh
set driver_script_dir = ${exp_setup_dir}

####################################################################
# Compile model (just once)
####################################################################
if ($compile_model > 0) then

   setenv CASE     $execase
   setenv CASEROOT ${PTMP}/cases/$CASE
   setenv RUNDIR   ${PTMP}/run/$CASE

   # Create and set up new case

   echo
   echo Start to create case
   echo
   source ${case_setup_script}
   echo
   echo Finished creating case
   echo

   # Build the model

   cd $CASEROOT
   ./xmlchange -file env_build.xml -id DEBUG   -val $debug

   echo
   echo Start to build model 
   echo
   if ( $postCIME <= 2 ) then
      ./$CASE.build
    else
      ./case.build
    endif

    echo 
    echo Finished building the model.
    echo

endif

#####################################################################
# Conduct simulation(s)
#####################################################################
if ($run_model > 0) then

  # initial condition files

  set atmInitList = (\
   ${initDir}'init_gen_clim_FC5_default.cam.i.2008-01-01-00000.nc' \
   ${initDir}'init_gen_clim_FC5_default.cam.i.2008-02-01-00000.nc' \
   ${initDir}'init_gen_clim_FC5_default.cam.i.2008-03-01-00000.nc' \
   ${initDir}'init_gen_clim_FC5_default.cam.i.2008-04-01-00000.nc' \
   ${initDir}'init_gen_clim_FC5_default.cam.i.2008-05-01-00000.nc' \
   ${initDir}'init_gen_clim_FC5_default.cam.i.2008-06-01-00000.nc' \
   ${initDir}'init_gen_clim_FC5_default.cam.i.2008-07-01-00000.nc' \
   ${initDir}'init_gen_clim_FC5_default.cam.i.2008-08-01-00000.nc' \
   ${initDir}'init_gen_clim_FC5_default.cam.i.2008-09-01-00000.nc' \
   ${initDir}'init_gen_clim_FC5_default.cam.i.2008-10-01-00000.nc' \
   ${initDir}'init_gen_clim_FC5_default.cam.i.2008-11-01-00000.nc' \
   ${initDir}'init_gen_clim_FC5_default.cam.i.2008-12-01-00000.nc' \
   )

   set ninit = $#atmInitList

   if ( $irE > $ninit ) then
      set irE = $ninit
   endif

   #---------------------------------
   cd $driver_script_dir

   set template = './cnvg_template_'$testconfig'.csh'

   #---------------------------------
   # GROUP and DTIME loops
   #---------------------------------
   set ig = $igS
   while ( $ig <= $igE )
    set group = ${groupList[$ig]}
    set nl_file = ${driver_script_dir}/"namelist_files/cam_nl_"${group}
    if ( ! -f $nl_file) then
     echo $nl_file " does not exist"
     exit 1
    endif

    set idtime = 1
    while ( $idtime <= $ndtime )

       set dtime  = $dtimeList[$idtime]
   
      #------------------------------------
      # create script for each realization
      #------------------------------------
      set ir = $irS
      while ( $ir <= $irE )
      
         set atm_init = $atmInitList[$ir]
         set irstring = `printf "%02d" ${ir}`
   
         set case_prefix     = ${group}_${RESOLUTION}_${COMPSET}_${COMPILER}_${MACH}
         set case_suffix     = "DT"`printf "%04d" ${dtime}`"_"${irstring}_${testconfig}
         set case = ${case_prefix}_${case_suffix}
         set caseroot = ${PTMP}/cases/${case_prefix}/${case_suffix}
         set rundir   = ${PTMP}/run/${case_prefix}/${case_suffix}
   
         set tmp_script = 'tmp_script_'`date +%F-%H%M%S-%N`
   
         cat <<EOF > $tmp_script
#!/bin/csh
date
setenv PROJECT         '$CESM_PROJ'
setenv CESM_PROJ       '$CESM_PROJ'
setenv CESM_EMAIL      '$CESM_EMAIL'

setenv MACH            '$MACH'
setenv COMPILER        '$COMPILER'
setenv RESOLUTION      '$RESOLUTION'
setenv NTASKS_PER_INST '$NTASKS_PER_INST'
setenv NINST           '$NINST'
setenv NTHRDS         '$NTHRDS'
         
setenv CCSMTAG  '$CCSMTAG'
setenv COMPSET  '$COMPSET'
setenv CCSMROOT '$CCSMROOT'
setenv PTMP     '$PTMP'
setenv EXEDIR   '$EXEDIR'
setenv CSMDATA  '$CSMDATA'
         
setenv CASE      '$case'
setenv CASEROOT  '$caseroot'
setenv RUNDIR    '$rundir'

setenv postCIME $postCIME

set case_setup_script = '$case_setup_script' 

set dtime    = $dtime
set atm_init = $atm_init

set ncycle   = 1       # 1 cycle in total, i.e., no restart

@ nlen = 3600 / $dtime   # run model for 60 minutes.
set stop_o = 'nsteps'

@ nhtfrq  = 1800 / $dtime    # output every 30 minutes
set mfilt = 1                # 1 time step per file 

set nl_file = $nl_file
EOF

      set jobscript = ${PTMP}/job_${case}.csh
      cat $tmp_script $template > $jobscript
      \rm $tmp_script

      echo Created heavy-wgt script$jobscript

      # run the heavy-wgt script to create case for a simulation
      csh $jobscript

      @ ir++
      end  #----------------------------------------------

   @ idtime++
   end  #----------------------------------------------

  @ ig++
  end  #----------------------------------------------

endif # if ($run_model > 1 )
