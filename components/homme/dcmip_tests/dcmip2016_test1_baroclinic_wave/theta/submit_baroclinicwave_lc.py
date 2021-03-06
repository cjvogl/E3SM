#!/usr/tce/bin/python

import os, sys

# Set MPI parameters based on machine
hostname = os.getenv('HOSTNAME')
if (hostname.startswith('cab')):
  num_nodes = 12
  num_procs_per_node = 16
  seconds_per_step = 0.1 # probably need to adjust this number
elif (hostname.startswith('quartz')):
  num_nodes = 5
  num_procs_per_node = 36
  seconds_per_step = 0.1
else:
  print('\n*** SCRIPT NOT IMPLEMENTED FOR THIS MACHINE **\n\n')
  sys.exit()

###############################################################

if (len(sys.argv) < 2):
  print('usage: ./submit_gravitywave.sh <executable name> [optional parameter name] [parameter value]')
  sys.exit()

# Set dictionary of default HOMME parameter names and values
paramDict = {
  'tsteptype':          '7',
  'tstep':            '320',
  'ndays':             '30',
  'ne':                 '8',
  'nu':            '3.0e16',
  'splitting':          '1',
  'rtol':          '1.0e-4',
  'atol':            '-1.0',
  'calcstats':        'true'
}

# Parse command line arguments and parameter dictionary
userList = []
for j in range(2,len(sys.argv),2):
  current = str(sys.argv[j])
  if (current in paramDict.keys()):
    paramDict[current] = sys.argv[j+1]
    userList.append(current)
  else:
    print('\n*** PARAMETER %s NOT VALID ***\n\n' % current)
    sys.exit()

# Use tstep_type if no user defined parameters given
if (not userList):
  userList = ['tsteptype']

# Set naming suffix
suffix = ''
for i,word in enumerate(userList):
  suffix += word + paramDict[word]
  if (i < len(userList)-1):
    suffix += '_'

# Create input namelist file text
namelist = \
"""&ctl_nl
nthreads               = 1
partmethod             = 4
topology               = "cube"
test_case              = "dcmip2016_test1"
ne                     = %s
qsize                  = 5
ndays                  = %s
statefreq              = 10
restartfreq            = -1
runtype                = 0
tstep                  = %s
integration            = "explicit"
tstep_type             = %s
rsplit                 = 1
qsplit                 = 1
nu                     = %s
nu_s                   = %s
nu_p                   = %s
nu_top                 = 0
limiter_option         = 8
hypervis_order         = 2
hypervis_subcycle      = 1
moisture               = 'dry'
theta_hydrostatic_mode = .false.
dcmip16_prec_type      = -1
dcmip16_pbl_type       = -1
/
&vert_nl
vform                  = "ccm"
vfile_mid              = "camm-30.ascii"
vfile_int              = "cami-30.ascii"
/
&analysis_nl
output_prefix          = "r400-dry-"
output_dir             = "./output_%s/"
output_timeunits       = 1
output_frequency       = 1
output_varnames1       = 'T','w'
interp_type            = 0
output_type            = 'netcdf'
num_io_procs           = 16
interp_nlon            = 360
interp_nlat            = 181
/
&arkode_nl
imex_splitting         = %s
rel_tol                = %s
abs_tol                = %s
calc_nonlinear_stats   = .%s.
/
&prof_inparm
profile_outpe_num      = 100
profile_single_file	   = .true.
/""" % (paramDict['ne'], paramDict['ndays'], paramDict['tstep'],
        paramDict['tsteptype'],paramDict['nu'],paramDict['nu'],paramDict['nu'],
        suffix,paramDict['splitting'],
        paramDict['rtol'], paramDict['atol'],paramDict['calcstats'])
os.system("echo '%s' > input_%s.nl" % (namelist, suffix))

# Create job script
steps_per_day = 60*60*24/int(paramDict['tstep'])
runtime = seconds_per_step*steps_per_day*float(paramDict['ndays'])
runtime = min(60*60*24,int(runtime))
hours = runtime/3600
runtime = runtime % 3600
minutes = runtime/60
seconds = runtime % 60
walltime = str(hours/10)+str(hours%10) + ':' + \
            str(minutes/10)+str(minutes%10) + ':' + \
            str(seconds/10)+str(seconds%10)

script = \
"""#!/bin/bash

#SBATCH -A climalg # Allocation
#SBATCH -p pbatch # pbatch or pdebug
#SBATCH -N %d # nodes
#SBATCH -t %s # Max walltime (hh:mm:ss)
#SBATCH -J gravitywave_%s # Job name
#SBATCH -o %s.out # stdout file
#SBATCH -e %s.err # stderr file

date
echo
pwd
echo

time srun -n %d %s < input_%s.nl
""" % (num_nodes,walltime,suffix,suffix,suffix,
       num_nodes*num_procs_per_node,sys.argv[1],suffix)
os.system("echo '%s' > run_%s.sh" % (script,suffix))

# Create output directory
os.mkdir('./output_%s' % suffix)

# Submit job script
os.system("sbatch ./run_%s.sh" % suffix)
