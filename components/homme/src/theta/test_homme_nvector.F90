module testing_mod
  implicit none
  ! dummy module so that elem(:) is possible (otherwise you need to explicitly
  ! define the shape of elem
contains

subroutine test_homme_nvector(elem,nets,nete,tl)

  use time_mod,           only: TimeLevel_t
  use HommeNVector,       only: NVec_t, MakeHommeNVector
  use element_mod,        only: element_t
  use parallel_mod,       only: parallel_t
  use physical_constants, only: dd_pi
  use dimensions_mod,     only: np, nlev
  use iso_c_binding

  implicit none
  integer, intent(in)                   :: nets, nete
  type (element_t), intent(in)          :: elem(:)
  type (TimeLevel_t), intent(in)        :: tl

  type(NVec_t), target          :: x, y
  type(c_ptr)                   :: x_C, y_C, z_C
  integer(C_INT)                :: ierr(2)
  real(C_DOUBLE)                :: cval

  !=======Internals ============

  print *, ""
  print *, "***************************"
  print *, "Running HOMME NVector Tests"
  print *, ""

  ! Test MakeHommeNVector (leave a registry space for ExtClone)
  print *, "Testing MakeHommeNVector"
  call MakeHommeNVector(elem, nets, nete, tl%np1, x, ierr(1))
  call MakeHommeNVector(elem, nets, nete, tl%n0, y, ierr(2))
  x_C = c_loc(x)
  y_C = c_loc(y)
  if (ierr(1) == 1 .or. ierr(2) == 1) then
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
  if (dabs(4.d0*dd_pi*6.d0*float(nlev) - cval) < 1.d-12) then
    print *, "success"
    print *, ""
  else
    print '("got ",f10.12," instead of  ",f10.12)', &
        cval, 4.d0*dd_pi*6.d0*float(nlev)
    print *, "test failed"
    return
  end if

  print *, "***************************"
  print *, ""

  ! Test FNVExtLinearSum
  print *, "Testing FNVExtLinearSum"
  call FNVExtLinearSum(1.d0,x_C,-1.d0,x_C,y_C)
  call FNVExtDotProd(y_C,y_C,cval)
  if (dabs(cval) < 1.d-12) then
    print *, "success"
    print *, ""
  else
    print '("got ",f10.12," instead of  ",f10.12)', &
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
  if (dabs(cval) < 1.d-12) then
    print *, "success"
    print *, ""
  else
    print '("got ",f10.12," instead of  ",f10.12)', &
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

end subroutine test_homme_nvector

end module testing_mod
