#!/usr/tce/bin/python

import os, sys

# Set MPI parameters based on machine
hostname = os.getenv('HOSTNAME')
if (hostname.startswith('cab')):
  num_nodes = 12
  num_procs_per_node = 16
  seconds_per_step = 3.0 # probably need to adjust this number
elif (hostname.startswith('quartz')):
  num_nodes = 5
  num_procs_per_node = 36
  seconds_per_step = 3.0
else:
  print('\n*** SCRIPT NOT IMPLEMENTED FOR HTIS MACHINE **\n\n')
  sys.exit()
###############################################################

if (len(sys.argv) < 2):
  print('usage: ./submit_gravitywave.sh <executable name> [optional parameter name] [parameter value]')
  sys.exit()

# Set dictionary of default HOMME parameter names and values
paramDict = {
  'tsteptype':          '5',
  'tstep':            '1.0',
  'nmax':            '3600',
  'ne':                '27',
  'nu':             '5.0e8',
  'splitting':          '1',
  'rtol':          '1.0e-8',
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
output_frequency = int(paramDict['nmax'])/10
namelist = \
"""&ctl_nl
nthreads               = 1
partmethod             = 4
topology               = "cube"
test_case              = "dcmip2012_test3"
ne                     = %s
qsize                  = 0
nmax                   = %s
statefreq              = 60
restartfreq            = -1
runtype                = 0
tstep                  = %s
rsplit                 = 0
tstep_type             = %s
integration            = "explicit"
theta_hydrostatic_mode = .false.
nu                     = %s
nu_p                   = %s
hypervis_order         = 2
hypervis_subcycle      = 1
rearth                 = 50969.76
omega                  = 0.0
/
&vert_nl
vform                  = "ccm"
vanalytic              = 1
vtop                   = 2.73919e-1
/
&analysis_nl
output_dir             = "./output_%s/"
output_timeunits       = 0
output_frequency       = %d
output_varnames1       = 'u', 'v', 'T'
interp_type            = 0
output_type            = 'netcdf'
num_io_procs           = 16
interp_nlat            = 128
interp_nlon            = 256
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
/""" % (paramDict['ne'], paramDict['nmax'], paramDict['tstep'],
        paramDict['tsteptype'],paramDict['nu'],paramDict['nu'],
        suffix,output_frequency,paramDict['splitting'],
        paramDict['rtol'], paramDict['atol'],paramDict['calcstats'])
os.system("echo '%s' > input_%s.nl" % (namelist, suffix))

# Create job script
runtime = seconds_per_step*float(paramDict['nmax'])
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
