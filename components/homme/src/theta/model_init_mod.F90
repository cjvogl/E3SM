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

end module
