!-----------------------------------------------------------------
! Daniel R. Reynolds
! SMU Mathematics
!
! This file contains all required subroutines for interfacing to
! the ARKode Fortran interface, as well as a small number of
! auxiliary routines to assist.  This uses a custom NVector
! interface, implemented in Fortran.
!
! The following subroutines are 'users' of the ARKode Fortran
! interface:
!    arkode_init -- initializes NVectors and ARKode for subsequent
!                   solves (only called once)
!
! The following user-supplied subroutines are required by the
! ARKode Fortran interface:
!    farkifun -- implements all implicit portions of the ODE
!                right-hand side function
!    farkefun -- implements all explicit portions of the ODE
!                right-hand side function
!
! The following user-supplied subroutines are optional within
! the ARKode Fortran interface:
!    farkjtimes -- implements a Jacobian-vector product, J*v,
!                  where J(u) is the Jacobian of farkifun
!                  with respect to u.
!    farkpset -- performs setup for the preconditioner matrix
!                (called infrequently)
!    farkpsol -- performs the preconditioner solve (called every
!                linear iteration)
!
! The following are 'helper' subroutines that I often use when
! interfacing to SUNDIALS solvers from Fortran:
!    farkdiags -- writes ARKode solver diagnostics to the screen
!
! Copyright 2017; all rights reserved
!=================================================================



subroutine arkode_init(elem, nets, nete, tl, y_C, ierr)
  !-----------------------------------------------------------------
  ! Description: arkode_init initializes the ARKode solver.
  !   Arguments:
  !     elem - (obj*, input) element objects
  !     nets - (int, input) starting index for elem array
  !     nete - (int, input) ending index for elem array
  !       tl - (obj, input) timelevel object
  !      y_C - (ptr(3), output) C pointers to NVec_t template solution vectors
  !     ierr - (int, output) return flag: 0=>success,
  !             1=>recoverable error, -1=>non-recoverable error
  !-----------------------------------------------------------------
  !======= Inclusions ===========
  use element_mod,      only: element_t
  use HommeNVector,     only: NVec_t, MakeHommeNVector
  use time_mod,         only: TimeLevel_t, tstep
  use iso_c_binding

  !======= Declarations =========
  implicit none

  ! calling variables
  type(element_t),   intent(in)  :: elem(nets:nete)
  type(TimeLevel_t), intent(in)  :: tl
  integer,           intent(in)  :: nets, nete
  type(c_ptr),       intent(out) :: y_C(3)
  integer(C_INT),    intent(out) :: ierr

  ! local variables
  type(NVec_t), target :: y(3)
  integer(C_INT)  :: idef, iatol, imex, precLR, gstype, maxl
  real*8          :: rpar(1), lintol, tstart, rtol, atol, rout(40)
  integer(C_LONG) :: lidef, ipar(1), iout(40)
  integer :: i


  ! Butcher table parameters for RK2
  integer :: S_RK2 = 2
  integer :: Q_RK2 = 3
  integer :: P_RK2 = 0
  real*8  :: C_RK2(2) = (/ 0.d0, 2.d0/3.d0 /)
  real*8  :: A_RK2(4) = (/ 0.d0, 0.d0, 2.d0/3.d0, 0.d0 /)
  real*8  :: B_RK2(2) = (/ 1.d0/4.d0, 3.d0/4.d0 /)

  !  -!!$       Ullrich 3rd order 5 stage
!  -!!$          s = 6;
!  -!!$          q = 3;
!  -!!$          A = zeros(s,s);
!  -!!$          A(2,1) = 0.2;
!  -!!$          A(3,2) = 0.2;
!  -!!$          A(4,3) = 1/3;
!  -!!$          A(5,4) = 2/3;
!  -!!$          A(6,1) = 0.25;  A(6,5) = 0.75;
!  -!!$          c = [0; 0.2; 0.2; 1/3; 2/3; 1];
!  -!!$          b = [0.25, 0, 0, 0, 0.75, 0];


  !======= Internals ============

  ! atol - absolute tolerance (for iterative solves)
  iatol = 1    ! specify type for atol: 1=scalar, 2=array
  atol = 1d-1    ! do all solution components have unit magnitude, or coul
               ! their units vary considerably?  If units can vary, then
               ! a scalar-valued atol is a **bad** idea

  ! rtol - relative tolerance (for iterative solves)
  rtol = 1d-1    ! do you really only want one digit of accuracy?  When us
               ! fixed time steps and an explicit method this input is u
               ! but in all other cases it corresponds to how tightly th
               ! are solved, and should rougly correspond with the desir
               ! fixed time steps and an explicit method this input is u
               ! but in all other cases it corresponds to how tightly th
               ! are solved, and should rougly correspond with the desir
               ! number of digits

  iout = 0
  rout = 0.d0
  tstart = 0.d0


  ! initialize error flag
  ierr = 0

  ! 'create' NVec_t objects that will correspond to the original 3 HOMME
  ! timelevels, assuming that tl%nm1, tl%n0, and tl%np1 are taken from the
  ! set {1,2,3}
  do i=1,3
    call MakeHommeNVector(elem, nets, nete, i, y(i), ierr)
    if (ierr /= 0) then
      print *,  'Error in MakeHommeNVector, ierr = ', ierr, '; halting'
      stop
    end if
    ! get C pointer
    y_C(i) = c_loc(y(i))
  end do

  ! initialize ARKode data & operators
  !    Nvector specs
  idef = 4    ! flag specifying which SUNDIALS solver will be used (4=ARKode)
  call fnvextinit(idef, ierr)
  if (ierr /= 0) then
     write(0,*) ' arkode_init: fnvextinit failed'
  endif

  !    ARKode dataspace
  imex = 1     ! specify problem type: 0=implicit, 1=explicit, 2=imex
  call farkmalloc(tstart, y_C(tl%n0), imex, iatol, rtol, atol, &
                  iout, rout, ipar, rpar, ierr)
  if (ierr /= 0) then
     write(0,*) ' arkode_init: farkmalloc failed'
  endif

  !      Indicate that we will set time step sizes ourselves, and the step
  !      size to use on the first time step (disable adaptivity)
  call farksetrin('FIXED_STEP', tstep, ierr)
  if (ierr /= 0) then
     write(0,*) ' arkode_init: farksetrin failed'
  endif


  !     Set ERK Butcher table
  call farkseterktable(S_RK2, Q_RK2, P_RK2, C_RK2, A_RK2, B_RK2, B_RK2, ierr)
  if (ierr /= 0) then
     write(0,*) ' arkode_init: farkseterktables failed'
  endif

  !      Requested order of accuracy for ARK method (3,4,5 are supported).
  !      Alternately, a user can supply a custom Butcher table pair to
  !      define their ARK method, by calling FARKSETARKTABLES()
  lidef = 4
  call farksetiin('ORDER', lidef, ierr)
  if (ierr /= 0) then
     write(0,*) ' arkode_init: farksetiin failed'
  endif


  !      To indicate that the implicit problem is linear, make the following
  !      call.  The argument specifies whether the linearly implicit problem
  !      changes as the problem evolves (1) or not (0)
!  lidef = 0
!  call farksetiin('LINEAR', lidef, ierr)
!  if (ierr /= 0) then
!     write(0,*) ' arkode_init: farksetiin failed'
!  endif

  !      Indicate use of the GMRES linear solver, the arguments indicate:
  !      precLR -- type of preconditioning: 0=none, 1=left, 2=right, 3=left+right
  !      gstype -- type of Gram-Schmidt orthogonalization: 1=modified, 2=classical
  !      maxl -- maximum size of Krylov subspace (# of iterations/vectors)
  !      lintol -- linear convergence tolerance factor (0 indicates default); this
  !                example is very stiff so it requires tight linear solves
!  precLR = 0
!  gstype = 1
!  maxl = 50
!  lintol = 1.d-3
!  call farkspgmr(precLR, gstype, maxl, lintol, ierr)
!  if (ierr /= 0) then
!     write(0,*) ' arkode_init: farkspgmr failed'
!  endif

  !      Indicate to use our own Jacobian-vector product routine (otherwise it
  !      uses a finite-difference approximation)
  !idef = 1
  !call farkspilssetjac(idef, ierr)
  !if (ierr /= 0) then
  !   write(0,*) ' arkode_init: farkspilssetjac failed'
  !endif

  !      Indicate to use our own preconditioner setup/solve routines (otherwise
  !      preconditioning is disabled)
  !idef = 1
  !call farkspilssetprec(idef, ierr)
  !if (ierr /= 0) then
  !   write(0,*) ' arkode_init: farkspilssetprec failed'
  !endif

  return
end subroutine arkode_init
!=================================================================



subroutine farkifun(t, y_C, fy_C, ipar, rpar, ierr)
  !-----------------------------------------------------------------
  ! Description: farkifun provides the implicit portion of the right
  !     hand side function for the ODE:   dy/dt = fi(t,y) + fe(t,y)
  !
  ! Arguments:
  !       t - (dbl, input) current time
  !     y_C - (ptr) C pointer to NVec_t containing current solution
  !    fy_C - (ptr) C pointer to NVec_t to hold right-hand side function
  !    ipar - (long int(*), input) integer user parameter data
  !           (passed back here, unused)
  !    rpar - (dbl(*), input) real user parameter data (passed here,
  !           unused)
  !    ierr - (int, output) return flag: 0=>success,
  !            1=>recoverable error, -1=>non-recoverable error
  !-----------------------------------------------------------------
  !======= Inclusions ===========
  use iso_c_binding
  use kinds,            only: real_kind
  use HommeNVector,     only: NVec_t
  use element_mod,      only: element_t
  use hybrid_mod,       only: hybrid_t
  use kinds,            only: real_kind
  use derivative_mod,   only: derivative_t
  use hybvcoord_mod,    only: hvcoord_t
  use dimensions_mod,   only: np,nlev
  use prim_advance_mod, only: compute_andor_apply_rhs, dt_save, eta_ave_w_save, &
                              qn0_save, hvcoord_ptr, hybrid_ptr, deriv_ptr

  !======= Declarations =========
  implicit none

  ! calling variables
  real*8,            intent(in)         :: t
  type(c_ptr),       intent(in), target :: y_C
  type(c_ptr),       intent(in), target :: fy_C
  integer(C_LONG),   intent(in)         :: ipar(1)
  real*8,            intent(in)         :: rpar(1)
  integer(C_INT),    intent(out)        :: ierr

  ! local variables
  type(NVec_t), pointer :: y => NULL()
  type(NVec_t), pointer :: fy => NULL()
  real (kind=real_kind) :: ci

  integer :: ie, inlev, inpx, inpy

  !======= Internals ============

  ! set return value to success
  ierr = 0

  ! dereference pointer for NVec_t objects
  call c_f_pointer(y_C, y)
  call c_f_pointer(fy_C, fy)

  ! determine 'stage time' from t, tcur and dt
  ! TODO: we don't have access to tcur here
  ! ci = (t - tcur)/dt_save
  ci = 1.d0 ! DEBUG

  ! set return value to success
  ierr = 0

  ! dereference pointer for NVec_t objects
  call c_f_pointer(y_C, y)
  call c_f_pointer(fy_C, fy)

  ! The function call to compute_andor_apply_rhs is as follows:
  !  compute_andor_apply_rhs(np1, nm1, n0, qn0, dt2, elem, hvcoord, hybrid, &
  !     deriv, nets, nete, compute_diagnostics, eta_ave_w, scale1, scale2, scale3)
  !
  !  This call returns the following:
  !
  !   u(np1) = scale3*u(nm1) + dt2*DSS[ nonstiffRHS(u(n0))*scale1 + stiffRHS(un0)*scale2 ]
  !
  !   nonstiffRHS and the stiffRHS are determined within the function and can be change by
  !   multiplying different terms by scale1 and scale2
  !
  !  Setting scale1=scale2=1.0, scale3=0.0, and dt2=1.0 returns the full rhs
  !
  !  DSS is the averaging procedure for the active and inactive nodes
  !


  call compute_andor_apply_rhs(fy%tl_idx, fy%tl_idx, y%tl_idx, qn0_save, &
       1.d0, y%elem, hvcoord_ptr, hybrid_ptr, deriv_ptr, y%nets, y%nete, &
       .false., ci*eta_ave_w_save, 0.d0, 1.d0, 0.d0)

  return
end subroutine farkifun
!=================================================================




subroutine farkefun(t, y_C, fy_C, ipar, rpar, ierr)
  !-----------------------------------------------------------------
  ! Description: farkefun provides the explicit portion of the right
  !     hand side function for the ODE:   dy/dt = fi(t,y) + fe(t,y)
  !
  ! Arguments:
  !       t - (dbl, input) current time
  !     y_C - (ptr) C pointer to NVec_t containing current solution
  !    fy_C - (ptr) C pointer to NVec_t to hold right-hand side function
  !    ipar - (long int(*), input) integer user parameter data
  !           (passed back here, unused)
  !    rpar - (dbl(*), input) real user parameter data (passed here,
  !           unused)
  !    ierr - (int, output) return flag: 0=>success,
  !            1=>recoverable error, -1=>non-recoverable error
  !-----------------------------------------------------------------
  use iso_c_binding
  use kinds,            only: real_kind
  use HommeNVector,     only: NVec_t
  use element_mod,      only: element_t
  use hybrid_mod,       only: hybrid_t
  use derivative_mod,   only: derivative_t
  use hybvcoord_mod,    only: hvcoord_t
  use dimensions_mod,   only: np,nlev
  use prim_advance_mod, only: compute_andor_apply_rhs, dt_save, eta_ave_w_save, &
                              qn0_save, hvcoord_ptr, hybrid_ptr, deriv_ptr

  !======= Declarations =========
  implicit none

  ! calling variables
  real*8,            intent(in)         :: t
  type(c_ptr),       intent(in), target :: y_C
  type(c_ptr),       intent(in), target :: fy_C
  integer(C_LONG),   intent(in)         :: ipar(1)
  real*8,            intent(in)         :: rpar(1)
  integer(C_INT),    intent(out)        :: ierr

  ! local variables
  type(NVec_t), pointer :: y => NULL()
  type(NVec_t), pointer :: fy => NULL()
  real (kind=real_kind) :: ci

  integer :: ie, inlev, inpx, inpy

  !======= Internals ============

  ! set return value to success
  ierr = 0

  ! dereference pointer for NVec_t objects
  call c_f_pointer(y_C, y)
  call c_f_pointer(fy_C, fy)

  ! determine 'stage time' from t, tcur and dt
  ! TODO: we don't have access to tcur here
  ! ci = (t - tcur)/dt_save
  ci = 1.d0 ! DEBUG

  ! The function call to compute_andor_apply_rhs is as follows:
  !  compute_andor_apply_rhs(np1, nm1, n0, qn0, dt2, elem, hvcoord, hybrid, &
  !     deriv, nets, nete, compute_diagnostics, eta_ave_w, scale1, scale2, scale3)
  !
  !  This call returns the following:
  !
  !   u(np1) = scale3*u(nm1) + dt2*DSS[ nonstiffRHS(u(n0))*scale1 + stiffRHS(un0)*scale2 ]
  !
  !   nonstiffRHS and the stiffRHS are determined within the function and can be change by
  !   multiplying different terms by scale1 and scale2
  !
  !  Setting scale1=scale2=1.0, scale3=0.0, and dt2=1.0 returns the full rhs
  !
  !  DSS is the averaging procedure for the active and inactive nodes
  !

  call compute_andor_apply_rhs(fy%tl_idx, fy%tl_idx, y%tl_idx, qn0_save, &
       1.d0, y%elem, hvcoord_ptr, hybrid_ptr, deriv_ptr, y%nets, y%nete, &
       .false., ci*eta_ave_w_save, 1.d0, 0.d0, 0.d0)

  return
end subroutine farkefun
!=================================================================



subroutine farkdiags(iout, rout)
  !-----------------------------------------------------------------
  ! Description: subroutine to output arkode diagnostics
  !
  ! Arguments:
  !        iout - (long int*, input) integer optional outputs
  !        rout - (dbl*, input) real optional outputs
  !-----------------------------------------------------------------
  !======= Inclusions ===========
  use iso_c_binding

  !======= Declarations =========
  implicit none

  integer(C_LONG), intent(in) :: iout(40)
  real*8,          intent(in) :: rout(40)

  !======= Internals ============

  ! general solver statistics
  print *, '   '
  print *,  ' ARKode output values:'
  print '(4x,A,i9)','Total internal steps taken =',iout(3)
  print '(4x,A,i9)','   stability-limited steps =',iout(4)
  print '(4x,A,i9)','    accuracy-limited steps =',iout(5)
  print '(4x,A,i9)','  internal steps attempted =',iout(6)
  print '(4x,A,i9)','Total explicit rhs calls   =',iout(7)
  print '(4x,A,i9)','Total implicit rhs calls   =',iout(8)

  print '(4x,A,i9)','Total nonlinear iterations =',iout(11)
  if (iout(13) > 0)  print '(4x,A,i9)','Total root function calls  =',iout(13)

  ! linear solver statistics
  print '(4x,A,i9)','Total linear iterations    =',iout(21)
  if (iout(9) > 0)  print '(4x,A,i9)','Num lin solver setup calls =',iout(9)
  if (iout(17) > 0) print '(4x,A,i9)','Num lin solver rhs calls   =',iout(17)
  if (iout(18) > 0) print '(4x,A,i9)','Num lin solver Jac calls   =',iout(18)
  if (iout(19) > 0) print '(4x,A,i9)','Num PSet routine calls     =',iout(19)
  if (iout(20) > 0) print '(4x,A,i9)','Num PSolve routine calls   =',iout(20)

  ! error statistics
  if (iout(10) > 0) print '(4x,A,i9)','Num error test failures    =',iout(10)
  if (iout(12) > 0) print '(4x,A,i9)','Num nonlin conv failures   =',iout(12)
  if (iout(22) > 0) print '(4x,A,i9)','Num linear conv failures   =',iout(22)


  ! general time-stepping information
  print '(4x,A,es12.5)','First internal step size   =',rout(1)
  print '(4x,A,es12.5)','Last internal step size    =',rout(2)
  print '(4x,A,es12.5)','Next internal step size    =',rout(3)
  print '(4x,A,es12.5)','Current internal time      =',rout(4)
  print '(4x,A,es12.5)','Suggested tol scale factor =',rout(5)
  print *, '   '

  return
end subroutine farkdiags
!=================================================================



!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
! NOTE: All of the remaining subroutines are not required when
! interfacing with ARKode, and so they are not currently
! implemented.  These may be provided to improve ARKode
! performance on this application; if so they should be 'enabled'
! by calling the relevant ARKode interface routine.
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


subroutine farkjtimes(v_C, Jv_C, t, y_C, fy_C, h, ipar, rpar, v1_C, ierr)
  !-----------------------------------------------------------------
  ! Description: farkjtimes provides the Jacobian-vector product
  !    routine for the linearized Newton system.
  !
  !  Arguments:
  !       v_C - (ptr) C pointer to NVec_t vector to multiply
  !      Jv_C - (ptr) C pointer to NVec_t for result of Jacobian-vector product
  !         t - (dbl, input) current time
  !       y_C - (ptr) C pointer to NVec_t containing current solution
  !      fy_C - (ptr) C pointer to NVec_t containing current implicit ODE rhs
  !         h - (dbl, input) time step size for last internal step
  !      ipar - (long int(*), input) integer user parameter data
  !             (passed back here, unused)
  !      rpar - (dbl(*), input) real user parameter data (passed here,
  !             unused)
  !      v1_C - (ptr) C pointer to NVec_t scratch vector
  !      ierr - (int, output) return flag: 0=>success,
  !             1=>recoverable error, -1=>non-recoverable error
  !-----------------------------------------------------------------
  !======= Inclusions ===========
  use iso_c_binding
  use HommeNVector,   only: NVec_t
  use dimensions_mod, only: np, nlev

  !======= Declarations =========
  implicit none

  ! calling variables
  type(c_ptr),     intent(in), target :: v_C
  type(c_ptr),     intent(in), target :: Jv_C
  real*8,          intent(in)         :: t
  type(c_ptr),     intent(in), target :: y_C
  type(c_ptr),     intent(in), target :: fy_C
  real*8,          intent(in)         :: h
  integer(C_LONG), intent(in)         :: ipar(1)
  real*8,          intent(in)         :: rpar(1)
  type(c_ptr),     intent(in), target :: v1_C
  integer(C_INT),  intent(out)        :: ierr

  ! local variables
  type(NVec_t), pointer :: v  => NULL()
  type(NVec_t), pointer :: Jv => NULL()
  type(NVec_t), pointer :: y  => NULL()
  type(NVec_t), pointer :: fy => NULL()
  type(NVec_t), pointer :: v1 => NULL()


  !======= Internals ============

  ! set return value to success
  ierr = 0

  ! dereference pointer for NVec_t objects
  call c_f_pointer(v_C, v)
  call c_f_pointer(Jv_C, Jv)
  call c_f_pointer(y_C, y)
  call c_f_pointer(fy_C, fy)
  call c_f_pointer(v1_C, v1)

  ! perform matrix-vector product, inserting result into Jv

  return
end subroutine farkjtimes
!=================================================================




subroutine farkpset(t, y_C, fy_C, jok, jcur, gamma, h, ipar, rpar, &
                    v1_C, v2_C, v3_C, ierr)
  !-----------------------------------------------------------------
  ! Description: farkpset provides the preconditioner setup
  !    routine for the linearized Newton system.
  !
  !  Arguments:
  !         t - (dbl, input) current time
  !       y_C - (ptr) C pointer to NVec_t containing current solution
  !      fy_C - (ptr) C pointer to NVec_t containing current implicit ODE rhs
  !       jok - (int, input) flag denoting whether to recompute
  !              Jacobian-related data: 0=>recompute, 1=>unnecessary
  !      jcur - (int, output) output flag to say if Jacobian data
  !              was recomputed: 1=>was recomputed, 0=>was not
  !     gamma - (dbl, input) the scalar appearing in the Newton matrix
  !              A = M-gamma*J
  !         h - (dbl, input) time step size for last internal step
  !      ipar - (long int(*), input) integer user parameter data
  !             (passed back here, unused)
  !      rpar - (dbl(*), input) real user parameter data (passed here,
  !             unused)
  !      v1_C - (ptr) C pointer to NVec_t scratch vector
  !      v2_C - (ptr) C pointer to NVec_t scratch vector
  !      v3_C - (ptr) C pointer to NVec_t scratch vector
  !      ierr - (int, output) return flag: 0=>success,
  !             1=>recoverable error, -1=>non-recoverable error
  !-----------------------------------------------------------------
  !======= Inclusions ===========
  use iso_c_binding
  use HommeNVector,   only: NVec_t
  use dimensions_mod, only: np, nlev

  !======= Declarations =========
  implicit none

  ! calling variables
  real*8,          intent(in)         :: t
  type(c_ptr),     intent(in), target :: y_C
  type(c_ptr),     intent(in), target :: fy_C
  integer(C_INT),  intent(in)         :: jok
  integer(C_INT),  intent(out)        :: jcur
  real*8,          intent(in)         :: gamma
  real*8,          intent(in)         :: h
  integer(C_LONG), intent(in)         :: ipar(1)
  real*8,          intent(in)         :: rpar(1)
  type(c_ptr),     intent(in), target :: v1_C
  type(c_ptr),     intent(in), target :: v2_C
  type(c_ptr),     intent(in), target :: v3_C
  integer(C_INT),  intent(out)        :: ierr

  ! local variables
  type(NVec_t), pointer :: y  => NULL()
  type(NVec_t), pointer :: fy => NULL()
  type(NVec_t), pointer :: v1 => NULL()
  type(NVec_t), pointer :: v2 => NULL()
  type(NVec_t), pointer :: v3 => NULL()

  !======= Internals ============

  ! initialize return value to success, jcur to not-recomputed
  ierr = 0
  jcur = 0

  ! dereference pointer for NVec_t objects
  call c_f_pointer(y_C, y)
  call c_f_pointer(fy_C, fy)
  call c_f_pointer(v1_C, v1)
  call c_f_pointer(v2_C, v2)
  call c_f_pointer(v3_C, v3)

  ! return if no preconditioner update is required
  if (jok == 1)  return

  ! update the preconditioner

  ! set Jacobian recomputation flag
  jcur = 1

  return
end subroutine farkpset
!=================================================================




subroutine farkpsol(t, y_C, fy_C, r_C, z_C, gamma, delta, lr, ipar, &
                    rpar, vt_C, ierr)
  !-----------------------------------------------------------------
  ! Description: farkpsol provides the preconditioner solve routine
  !    for the preconditioning of the linearized Newton system.
  !         i.e. solves P*z = r
  !
  ! Arguments:
  !         t - (dbl, input) current time
  !       y_C - (ptr) C pointer to NVec_t containing current solution
  !      fy_C - (ptr) C pointer to NVec_t containing current implicit ODE rhs
  !       r_C - (ptr) C pointer to NVec_t rhs vector of prec. system
  !       z_C - (ptr) C pointer to NVec_t solution vector of prec. system
  !     gamma - (dbl, input) scalar appearing in the Newton Matrix
  !             A = M-gamma*J
  !     delta - (dbl, input) desired tolerance if using an iterative
  !             method.  In that case, solve until
  !                  Sqrt[Sum((r-Pz).*ewt)^2] < delta
  !             where the ewt vector is obtainable by calling
  !             FARKGETERRWEIGHTS()
  !        lr - (int, input) flag indicating preconditioning type to
  !             apply: 1 => left,  2 => right
  !      ipar - (long int(*), input) integer user parameter data
  !             (passed back here, unused)
  !      rpar - (dbl(*), input) real user parameter data (passed here,
  !             unused)
  !      vt_C - (ptr) C pointer to NVec_t scratch vector
  !      ierr - (int, output) return flag: 0=>success,
  !             1=>recoverable error, -1=>non-recoverable error
  !-----------------------------------------------------------------
  !======= Inclusions ===========
  use iso_c_binding
  use HommeNVector,   only: NVec_t
  use dimensions_mod, only: np, nlev

  !======= Declarations =========
  implicit none

  ! calling variables
  real*8,          intent(in)         :: t
  type(c_ptr),     intent(in), target :: y_C
  type(c_ptr),     intent(in), target :: fy_C
  type(c_ptr),     intent(in), target :: r_C
  type(c_ptr),     intent(in), target :: z_C
  real*8,          intent(in)         :: gamma
  real*8,          intent(in)         :: delta
  integer(C_INT),  intent(in)         :: lr
  integer(C_LONG), intent(in)         :: ipar(1)
  real*8,          intent(in)         :: rpar(1)
  type(c_ptr),     intent(in), target :: vt_C
  integer(C_INT),  intent(out)        :: ierr

  ! local variables
  type(NVec_t), pointer :: y  => NULL()
  type(NVec_t), pointer :: fy => NULL()
  type(NVec_t), pointer :: r  => NULL()
  type(NVec_t), pointer :: z  => NULL()
  type(NVec_t), pointer :: vt => NULL()

  !======= Internals ============

  ! set return value to success
  ierr = 0

  ! dereference pointer for NVec_t objects
  call c_f_pointer(y_C, y)
  call c_f_pointer(fy_C, fy)
  call c_f_pointer(r_C, r)
  call c_f_pointer(z_C, z)
  call c_f_pointer(vt_C, vt)

  ! perform preconditioner solve to fill z

  return
end subroutine farkpsol
!=================================================================
