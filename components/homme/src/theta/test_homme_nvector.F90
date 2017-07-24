subroutine test_homme_nvector(elem,nets,nete,tl)

  use time_mod,         only: TimeLevel_t
  use HommeNVector,     only: NVec_t, MakeHommeNVector
  use element_mod,      only: element_t
  use parallel_mod,     only: parallel_t
  use dimensions_mod, only: np, nlev
  use, intrinsic :: iso_c_binding

  implicit none
  type (element_t), intent(in), target  :: elem(:)
  type (TimeLevel_t), intent(in)        :: tl
  integer, intent(in)                   :: nets, nete

  type(NVec_t), target          :: x
  type(c_ptr)                   :: x_C
  integer(C_INT)                :: ierr

  !=======Internals ============

  print *, ""
  print *, "***************************"
  print *, "Running HOMME NVector Tests"
  print *, ""

  ! Test MakeHommeNVector
  print *, "Testing MakeHommeNVector"
  call MakeHommeNVector(elem, nets, nete, tl%n0, x, ierr)
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
