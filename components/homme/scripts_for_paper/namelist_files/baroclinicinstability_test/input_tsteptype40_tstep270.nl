&ctl_nl
theta_hydrostatic_mode = .false.
dcmip4_moist           = 0
dcmip4_X               = 1.0
NThreads               = 1
partmethod             = 4
topology               = "cube"
test_case              = "dcmip2012_test4"
u_perturb              = 1
rotate_grid            = 0
ne                     = 30
qsize                  = 0
nmax                   = 4800
statefreq              = 100
runtype                = 0
mesh_file              = "/dev/null"
tstep                  = 270
rsplit                 = 1
qsplit                 = 1
tstep_type             = 40
integration            = "explicit"
nu                     = 1.0e15
nu_div                 = 1.0e15
nu_p                   = 1.0e15
nu_q                   = 1.0e15
nu_s                   = 1.0e15
nu_top                 = 0
se_ftype               = 0
limiter_option         = 9
vert_remap_q_alg       = 0
hypervis_scaling       = 0
hypervis_order         = 2
hypervis_subcycle      = 3
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
output_timeunits       = 0
output_frequency       = 320
output_start_time      = 0
output_varnames1       = u,v,w,T,geo
num_io_procs           = 16
output_type            = netcdf
output_dir             = "./output_tsteptype40_tstep270/"
/
&arkode_nl
rel_tol                = 1.0e-6
abs_tol                = -1
calc_nonlinear_stats   = .true.
use_column_solver      = .true.
/
