#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

!
!  model specific initialization code called from prim_init2
!  (threaded initialization code)
!
!  most models do nothing.  introduced for preqx_acc to initialize
!  GPU related data
!
module model_init_mod
  implicit none

contains

  subroutine model_init2( elem , deriv, nets, nete, tl )
    use element_mod, only: element_t
    use derivative_mod, only: derivative_t
    use time_mod,       only: TimeLevel_t, tstep
    use prim_advance_mod, only: y_C, ynew_C
    use HommeNVector, only: NVec_t, MakeHommeNVector
    use iso_c_binding
    implicit none
    type(element_t),          intent(in) :: elem(:)
    type(derivative_t),       intent(in) :: deriv
    type(TimeLevel_t),        intent(in) :: tl
    integer,                  intent(in) :: nets
    integer,                  intent(in) :: nete

    type(NVec_t), target :: y, ynew
    real*8 :: tstart, atol, rtol, rout(40)
    integer(C_LONG) :: iout(40)
    integer(C_INT) :: ierr

#ifdef TEST_HOMME_NVEC_INLINE
    call test_homme_nvector(elem,nets,nete,tl)
    stop
#endif

    ! 'create' NVec_t objects 'y' and 'ynew' to hold current solution
    call MakeHommeNVector(elem, nets, nete, tl%n0, y, ierr)
    if (ierr /= 0) then
      print *,  'Error in MakeHommeNVector, ierr = ', ierr, '; halting'
      stop
    end if
    call MakeHommeNVector(elem, nets, nete, tl%np1, ynew, ierr)
    if (ierr /= 0) then
      print *,  'Error in MakeHommeNVector, ierr = ', ierr, '; halting'
      stop
    end if

    ! create C pointers to y and ynew
    y_C = c_loc(y)
    ynew_C = c_loc(ynew)

    ! initialize arkode
    atol = 1d-1    ! do all solution components have unit magnitude, or coul
                   ! their units vary considerably?  If units can vary, then
                   ! a scalar-valued atol is a **bad** idea
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
     call arkode_init(tstart, tstep, y_C, rtol, atol, iout, rout, ierr)
     if (ierr /= 0) then
       print *,  'Error in arkode_init, ierr = ', ierr, '; halting'
       stop
     end if
       !! Set additional ARKode options (e.g. Butcher table)
       !!$       CALL FARKSETERKTABLE(S, Q, P, C, A, B, B2, IER)
       !!$
       !!$     The arguments are:
          !  -!!$       S = the number of stages in the table [int, input]
          !  -!!$       Q = the global order of accuracy of the method [int, input]
          !  -!!$       P = the global order of accuracy of the embedding [int, input]
          !  -!!$       C = array of length S containing the stage times [realtype, input]
          !  -!!$       A = array of length S*S containing the ERK coefficients (stored in
          !  -!!$           row-major, "C", order) [realtype, input]
          !  -!!$       B = array of length S containing the solution coefficients
          !  -!!$           [realtype, input]
          !  -!!$       B2 = array of length S containing the embedding coefficients
          !  -!!$           [realtype, input]
          !  -
          !  -!!$       RK2:
          !  -!!$          A = [ 0, 0; 2/3, 0];   !!! A = (/ 0.d0, 0.d0, 2.d0/3.d0, 0.d0 /)
          !  -!!$          b = [ 1/4, 3/4];
          !  -!!$          c = [ 0; 2/3];
          !  -!!$          q = 3;  p = 0;
          !  -!!$          b2 = [ 1/4, 3/4 ];
          !  -!!$          s = 2;
          !  -
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
          !  -
          !  +


  end subroutine

#ifdef TEST_HOMME_NVEC_INLINE
  subroutine test_homme_nvector(elem,nets,nete,tl)

    use time_mod,           only: TimeLevel_t
    use HommeNVector,       only: NVec_t, MakeHommeNVector, SetHommeNVectorComm
    use element_mod,        only: element_t
    use parallel_mod,       only: parallel_t
    use physical_constants, only: dd_pi
    use dimensions_mod,     only: np, nlev
    use iso_c_binding

    implicit none
#include <mpif.h>
    integer, intent(in)                   :: nets, nete
    type (element_t), intent(in)          :: elem(:)
    type (TimeLevel_t), intent(in)        :: tl

    type(NVec_t), target          :: x, y, z
    type(c_ptr)                   :: x_C, y_C, z_C
    integer(C_INT)                :: ierr1, ierr2
    real(C_DOUBLE)                :: cval, tol, cval2

    !=======Internals ============

    tol = 1.d-12

    call SetHommeNVectorComm(MPI_COMM_WORLD)

    print *, ""
    print *, "***************************"
    print *, "Running HOMME NVector Tests"
    print *, ""

    ! Test MakeHommeNVector (leave a registry space for ExtClone)
    print *, "Testing MakeHommeNVector"
    call MakeHommeNVector(elem, nets, nete, tl%np1, x, ierr1)
    call MakeHommeNVector(elem, nets, nete, tl%n0, y, ierr2)
    x_C = c_loc(x)
    y_C = c_loc(y)
    if (ierr1 == 1 .or. ierr2 == 1) then
      stop "Error making Homme NVector"
    else
      print *, "success"
      print *, ""
    end if

    ! Test FNVExtPrint (output should have some initial condition in it)
    print *, "Testing FNVExtPrint"
    call FNVExtPrint(y_C)
    print *, "if values contain some initial condition data, success"
    print *, ""

    print *, "***************************"
    print *, ""

    ! Test FNVExtConst
    print *, "Testing FNVExtConst"
    call FNVExtConst(1.d0,y_C)
    call FNVExtPrint(y_C)
    print *, "if all values are 1, success"
    print *, ""

    print *, "***************************"
    print *, ""

    ! Test FNVExtDotProd
    print *, "Testing FNVExtDotProd"
    call FNVExtDotProd(y_C,y_C,cval)
    if (dabs(4.d0*dd_pi*6.d0*nlev - cval) < tol) then
      print *, "success"
      print *, ""
    else
      print '(" got ",d20.12," instead of  ",d20.12)', &
          cval, 4.d0*dd_pi*6.d0*nlev
      print *, "test failed"
      return
    end if

    print *, "***************************"
    print *, ""

    ! Test FNVExtLinearSum
    print *, "Testing FNVExtLinearSum"
    call FNVExtLinearSum(1.d0,x_C,-1.d0,x_C,y_C)
    call FNVExtDotProd(y_C,y_C,cval)
    if (dabs(cval) < tol) then
      print *, "success"
      print *, ""
    else
      print '(" got ",d20.12," instead of  ",d20.12)', &
          cval, 0.d0
      print *, "test failed"
      return
    end if

    print *, "***************************"
    print *, ""

    ! Test FNVExtClone
    print *, "Testing FNVExtClone"
    call FNVExtClone(x_C,z_C)
    if (.not. c_associated(z_C)) then
      print *, "test failed (FNVExtClone returned C_NULL_PTR)"
      return
    end if
    call FNVExtLinearSum(1.d0,x_C,-1.d0,z_C,y_C)
    call FNVExtDotProd(y_C,y_C,cval)
    if (dabs(cval) < tol) then
      print *, "success"
      print *, ""
    else
      print '(" got ",d20.12," instead of  ",d20.12)', &
          cval, 0.d0
      print *, "test failed (FNVExtDotProd returned non-zero value)"
      return
    end if

    print *, "***************************"
    print *, ""

    ! Test FNVExtDestroy
    print *, "Testing FNVExtDestroy"
    call FNVExtDestroy(z_C)
    print *, "success"

    print *, "***************************"
    print *, ""

    ! Test FNVExtProd
    print *, "Testing FNVExtProd"
    call FNVExtConst(2.d0,y_C)
    call FNVExtProd(x_C,y_C,y_C)
    call FNVExtLinearSum(2.d0,x_C,-1.d0,y_C,y_C)
    call FNVExtDotProd(y_C,y_C,cval)
    if (dabs(cval) < tol) then
      print *, "success"
      print *, ""
    else
      print '(" got ",d20.12," instead of  ",d20.12)', &
          cval, 0.d0
      print *, "test failed"
      return
    end if

    print *, "***************************"
    print *, ""

    ! Test FNVAddConst
    print *, "Testing FNVAddConst"
    call FNVExtAddConst(5.d0,x_C,y_C)
    call FNVExtAddConst(-3.d0,y_C,y_C)
    call FNVExtAddConst(-2.d0,y_C,y_C)
    call FNVExtLinearSum(1.d0,x_C,-1.d0,y_C,y_C)
    call FNVExtDotProd(y_C,y_C,cval)
    if (dabs(cval) < tol) then
      print *, "success"
      print *, ""
    else
      print '(" got ",d20.12," instead of  ",d20.12)', &
          cval, 0.d0
      print *, "test failed"
      return
    end if

    print *, "***************************"
    print *, ""


    ! Test FNVExtDiv
    print *, "Testing FNVExtDiv"
    call FNVExtAddConst(tol,x_C,y_C)
    call FNVExtDiv(x_C,y_C,y_C)
    call FNVExtProd(x_C,y_C,y_C)
    call FNVExtLinearSum(-1.d0,x_C,1.d0,y_C,y_C)
    call FNVExtDotProd(y_C,y_C,cval)
    if (dabs(cval) < tol) then
      print *, "success"
      print *, ""
    else
      print '(" got ",d20.12," instead of  ",d20.12)', &
          cval, 0.d0
      print *, "test failed"
      return
    end if

    print *, "***************************"
    print *, ""

    ! Test FNVExtProd
    print *, "Testing FNVExtScale"
    call FNVExtScale(3.d0,x_C,y_C)
    call FNVExtLinearSum(-1.d0,y_C,3.d0,x_C,y_C)
    call FNVExtDotProd(y_C,y_C,cval)
    if (dabs(cval) < tol) then
      print *, "success"
      print *, ""
    else
      print '(" got ",d20.12," instead of  ",d20.12)', &
             cval, 0.d0
      print *, "test failed"
      return
    end if

    print *, "***************************"
    print *, ""

    ! Test FNVExtAbs
    print *, "Testing FNVExtAbs"
    call FNVExtConst(-2.d0,y_C)
    call FNVExtAbs(y_C,y_C)
    call FNVExtLinearSum(1.d0,y_C,-1.d0,y_C,y_C)
    call FNVExtDotProd(y_C,y_C,cval)
    if (dabs(cval) < tol) then
      print *, "success"
      print *, ""
    else
      print '(" got ",d20.12," instead of  ",d20.12)', &
             cval, 0.d0
      print *, "test failed"
      return
    end if

    print *, "***************************"
    print *, ""

    ! Test FNVExtInv
    print *, "Testing FNVExtInv"
    call FNVExtConst(4.d0,y_C)
    call FNVExtInv(y_C,y_C)
    call MakeHommeNVector(elem,nets,nete,tl%nm1,z,ierr1)
    z_C = c_loc(z)
    call FNVExtConst(0.25d0,z_C)
    call FNVExtLinearSum(1.d0,y_C,-1.d0,z_C,y_C)
    call FNVExtDotProd(y_C,y_C,cval)
    if (dabs(cval) < tol) then
      print *, "success"
      print *, ""
    else
      print '(" got ",d20.12," instead of  ",d20.12)', &
             cval, 0.d0
      print *, "test failed"
      return
    end if

    print *, "***************************"
    print *, ""

    ! Test FNVExtMaxNorm
    print *, "Testing FNVExtMaxNorm"
    call FNVExtMaxNorm(x_C,cval)
    x%elem((nets+nete)/2)%state%w(1,1,1,x%tl_idx) = cval + 1.d0
    call FNVExtMaxNorm(x_C,cval2)
    if (dabs(cval2-cval-1.d0) < tol) then
      print *, "success"
      print *, ""
    else
      print '(" got ",d20.12," instead of  ",d20.12)', &
             cval2, cval+1.d0
      print *, "test failed"
      return
    end if

    print *, "***************************"
    print *, ""

    ! Test FNVExtMin
    print *, "Testing FNVExtMin"
    call FNVExtConst(1.d0,y_C)
    y%elem((nets+nete)/2)%state%w(1,1,1,y%tl_idx) = 0.5d0
    call FNVExtMin(y_C,cval)
    if (dabs(cval-0.5d0) < tol) then
      print *, "success"
      print *, ""
    else
      print '(" got ",d20.12," instead of  ",d20.12)', &
             cval, 0.5d0
      print *, "test failed"
      return
    end if

    print *, "***************************"
    print *, ""

    ! Test FNVExtWrmsNorm
    print *, "Testing FNVExtWrmsNorm"
    call FNVExtConst(5.d0,y_C)
    call FNVExtConst(0.2d0,z_C)
    call FNVExtWrmsNorm(y_C,z_C,cval)
    if (dabs(cval-1.d0) < tol) then
      print *, "success"
      print *, ""
    else
      print '(" got ",d20.12," instead of  ",d20.12)', &
             cval, 1.0d0
      print *, "test failed"
      return
    end if

    print *, "***************************"
    print *, ""

  end subroutine test_homme_nvector
#endif

end module
