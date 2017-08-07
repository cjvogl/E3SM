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
    use prim_advance_mod, only: y_C, ynew_C
    use control_mod, only: tstep_type
    use element_mod, only: element_t
    use derivative_mod, only: derivative_t
    use time_mod,       only: TimeLevel_t, tstep
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

    if (tstep_type == 8) then
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
     end if

  end subroutine

end module
