subroutine test_homme_nvector(elem,hvcoord,hybrid,nets,nete,tl,par)

  ! It is assumed that deriv1 in prim_driver_base has been set before the call
  ! to test_homme_nvector().  Options are through 
  ! /src/share/prim_driver_base%prim_init1() or /src/share/prim_driver_base%prim_init2
  use prim_driver_base, only: deriv1
  use control_mod,      only: qsplit ! constant defined in /src/share/control_mod.F90
  use time_mod,         only: TimeLevel_t, TimeLevel_Qdp
  use HommeNVector,     only: NVec_t, MakeHommeNVector
  use element_mod,      only: element_t
  use hybvcoord_mod,    only: hvcoord_t
  use hybrid_mod,       only: hybrid_t
  use derivative_mod,   only: derivative_t
  use parallel_mod,     only: parallel_t
  use dimensions_mod, only: np, nlev
  use, intrinsic :: iso_c_binding

  implicit none
  type (element_t), intent(in), target  :: elem(:)
  type (hvcoord_t), intent(in)          :: hvcoord
  type (hybrid_t), intent(in)           :: hybrid
  type (TimeLevel_t), intent(in)        :: tl 
  type (parallel_t), intent(in)         :: par
  integer, intent(in)                   :: nets, nete

  type(NVec_t), target          :: x
  type(c_ptr)                   :: x_C
  integer                       :: ier, qn0, tl_idx

  !=======Internals ============

  print *, ""
  print *, "***************************"
  print *, "Running HOMME NVector Tests"
  print *, ""

  tl_idx = tl%n0
  call TimeLevel_Qdp(tl, qsplit, qn0)

  ! Test MakeHommeNVector
  print *, "Testing MakeHommeNVector"
  call MakeHommeNVector(elem, hvcoord, hybrid, deriv1, nets, nete, qn0, tl_idx, par, x, ier)  
  x_C = c_loc(x)
  if (ier == 1) then
    stop "Error making Homme NVector for x"
  else
    print *, "success"
    print *, ""
  end if

  ! Test FNVExtPrint
  print *, "Testing FNExtPrint"
  call FNVExtPrint(x_C)

  print *, "***************************"
  print *, ""

end subroutine test_homme_nvector
