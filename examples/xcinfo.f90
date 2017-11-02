program xcinfo
  use xc_f90_types_m
  use xc_f90_lib_m

  implicit none

  TYPE(xc_f90_pointer_t) :: xc_func
  TYPE(xc_f90_pointer_t) :: xc_info
  integer :: i
  character(len=120) :: s1, s2

  call xc_f90_func_init(xc_func, xc_info, XC_GGA_C_PBE, XC_UNPOLARIZED)

  select case(xc_f90_info_kind(xc_info))
  case(XC_EXCHANGE)
    write(*, '(a)') 'Exchange'
  case(XC_CORRELATION)
    write(*, '(a)') 'Correlation'
  case(XC_EXCHANGE_CORRELATION)
    write(*, '(a)') 'Exchange-correlation'
  case(XC_KINETIC)
    write(*, '(a)') 'Kinetic'
  end select

  call xc_f90_info_name(xc_info, s1)
  select case(xc_f90_info_family(xc_info))
  case (XC_FAMILY_LDA);       write(s2,'(a)') "LDA"
  case (XC_FAMILY_GGA);       write(s2,'(a)') "GGA"
  case (XC_FAMILY_HYB_GGA);   write(s2,'(a)') "Hybrid GGA"
  case (XC_FAMILY_MGGA);      write(s2,'(a)') "MGGA"
  case (XC_FAMILY_HYB_MGGA);  write(s2,'(a)') "Hybrid MGGA"
  case (XC_FAMILY_LCA);       write(s2,'(a)') "LCA"
  end select
  write(*, '(4a)') trim(s1), ' (', trim(s2), ')'

  i = 0
  call xc_f90_info_refs(xc_info, i, s1)
  do while(i >= 0)
    write(*, '(a,i1,2a)') '[', i, '] ', trim(s1)
    call xc_f90_info_refs(xc_info, i, s1)
  end do

  call xc_f90_func_end(xc_func)

end program xcinfo
