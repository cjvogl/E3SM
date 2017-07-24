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

  subroutine model_init2( elem , deriv, nets, nete, tl)
    use element_mod, only: element_t
    use derivative_mod, only: derivative_t
    use time_mod,       only: TimeLevel_t
    type(element_t)   , intent(in) :: elem(:)
    type(derivative_t), intent(in) :: deriv
    type(TimeLevel_t),  intent(in) :: tl
    integer,            intent(in) :: nets
    integer,            intent(in) :: nete

    ! do nothing
  end subroutine

end module
