&ctl_nl
nthreads               = 1
partmethod             = 4
topology               = "cube"
test_case              = "dcmip2012_test3"
ne                     = 27
qsize                  = 0
nmax                   = 600
statefreq              = 100
restartfreq            = -1
runtype                = 0
tstep                  = 0.5
rsplit                 = 0
tstep_type             = 42
integration            = "explicit"
theta_hydrostatic_mode = .false.
nu                     = 0.0
nu_p                   = 0.0
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
output_dir             = "./output_tsteptype42_tstep0.5_nmax600_nu0.0/"
output_timeunits       = 0
output_frequency       = 60
output_varnames1       = u,v,w,T,geo
interp_type            = 0
output_type            = netcdf
num_io_procs           = 16
interp_nlat            = 128
interp_nlon            = 256
/
&arkode_nl
rel_tol                = 1.0e-6
abs_tol                = -1.0
calc_nonlinear_stats   = .true.
use_column_solver      = .true.
/
&prof_inparm
profile_outpe_num      = 100
profile_single_file	   = .true.
/
