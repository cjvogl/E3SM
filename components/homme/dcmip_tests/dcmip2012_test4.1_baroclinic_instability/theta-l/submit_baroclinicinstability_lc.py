#!/usr/tce/bin/python

import os, sys

# Set MPI parameters based on machine
hostname = os.getenv('HOSTNAME')
if (hostname.startswith('cab')):
  num_nodes = 12
  num_procs_per_node = 16
  seconds_per_step = 1.0 # probably need to adjust this number
elif (hostname.startswith('quartz')):
  num_nodes = 5
  num_procs_per_node = 36
  seconds_per_step = 1.0
else:
  print('\n*** SCRIPT NOT IMPLEMENTED FOR THIS MACHINE **\n\n')
  sys.exit()

###############################################################

if (len(sys.argv) < 2):
  print('usage: ./submit_baroclinicinstability.sh <executable name> [optional parameter name] [parameter value]')
  sys.exit()

# Set dictionary of default HOMME parameter names and values
paramDict = {
  'tsteptype':          '7',
  'tstep':            '120',
  'ndays':             '30',
  'rsplit':            '10',
  'ne':                '20',
  'nu':            '3.7e15',
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
theta_hydrostatic_mode = .false.
dcmip4_moist           = 0
dcmip4_X               = 1.0
NThreads               = 1
partmethod             = 4
topology               = "cube"
test_case              = "dcmip2012_test4"
u_perturb              = 1
rotate_grid            = 0
ne                     = %s
qsize                  = 0
ndays                  = %s
statefreq              = 60
runtype                = 0
mesh_file              = "/dev/null"
tstep                  = %s
rsplit                 = %s
qsplit                 = 1
tstep_type             = %s
integration            = "explicit"
nu                     = %s
nu_div                 = %s
nu_p                   = %s
nu_q                   = %s
nu_s                   = %s
nu_top                 = 0
se_ftype               = 0
limiter_option         = 9
vert_remap_q_alg       = 0
hypervis_scaling       = 0
hypervis_order         = 2
hypervis_subcycle      = 1
/
&vert_nl
vform                  = "ccm"
vfile_mid              = "camm-30.ascii"
vfile_int              = "cami-30.ascii"
/
&prof_inparm
profile_outpe_num      = 100
profile_single_file	   = .true.
/
&analysis_nl
interp_gridtype        = 2
output_timeunits       = 1
output_frequency       = 1
output_start_time      = 7
output_varnames1       = 'ps','zeta','u','v','T'
num_io_procs           = 16
output_type            = 'netcdf'
output_prefix          = "nonhydro-X1-"
output_dir             = "./output_%s/"
/
&arkode_nl
imex_splitting         = %s
rel_tol                = %s
abs_tol                = %s
calc_nonlinear_stats   = .%s.
/""" % (paramDict['ne'], paramDict['ndays'], paramDict['tstep'], paramDict['rsplit'],
        paramDict['tsteptype'], paramDict['nu'], paramDict['nu'], paramDict['nu'], 
        paramDict['nu'], paramDict['nu'], suffix,paramDict['splitting'],
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
#SBATCH -J baroclinicinstability_%s # Job name
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
