program lxctest
  use xc_f90_types_m
  use xc_f90_lib_m

  implicit none

  TYPE(xc_f90_pointer_t) :: xc_func
  TYPE(xc_f90_pointer_t) :: xc_info
  real(8) :: rho(5) = (/0.1, 0.2, 0.3, 0.4, 0.5/)
  real(8) :: sigma(5) = (/0.2, 0.3, 0.4, 0.5, 0.6/)
  real(8) :: exc(5)
  integer :: i, vmajor, vminor, vmicro, func_id = 1

  call xc_f90_version(vmajor, vminor, vmicro)
  write(*,'("Libxc version: ",I1,".",I1,".",I1)') vmajor, vminor, vmicro

  call xc_f90_func_init(xc_func, xc_info, func_id, XC_UNPOLARIZED)

  select case (xc_f90_info_family(xc_info))
  case(XC_FAMILY_LDA)
    call xc_f90_lda_exc(xc_func, 5, rho(1), exc(1))
  case(XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
    call xc_f90_gga_exc(xc_func, 5, rho(1), sigma(1), exc(1))
  end select

  do i = 1, 5
    write(*,"(F8.6,1X,F9.6)") rho(i), exc(i)
  end do

  call xc_f90_func_end(xc_func)

end program lxctest
