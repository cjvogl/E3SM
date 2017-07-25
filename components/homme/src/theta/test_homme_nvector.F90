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

  type(NVec_t), target          :: x, y, z
  type(c_ptr)                   :: x_C, y_C, z_C
  integer(C_INT)                :: ierr1, ierr2
  real(C_DOUBLE)                :: cval, tol, cval2

  !=======Internals ============

  tol = 1.d-12

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

end module testing_mod
