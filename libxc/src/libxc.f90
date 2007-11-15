!! Copyright (C) 2003-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
!!
!! This program is free software; you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation; either version 2, or (at your option)
!! any later version.
!!
!! This program is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with this program; if not, write to the Free Software
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.
!!
!! $Id$

module xc_types_m

  type xc_func_t
    private
    integer, pointer :: func
  end type xc_func_t

  type xc_info_t
    private
    integer, pointer :: info
  end type xc_info_t

end module xc_types_m

module libxc_m
  
  use xc_types_m
  use libxc_funcs_m

  implicit none

  public
  public ::                             &
    xc_func_t,                          &
    xc_info_t,                          &
    xc_f90_info_number,                 &
    xc_f90_info_kind,                   &
    xc_f90_info_name,                   &
    xc_f90_info_family,                 &
    xc_f90_info_refs,                   &
    xc_f90_family_from_id,              &
    xc_f90_lda_init,                    &
    xc_f90_lda,                         &
    xc_f90_lda_vxc,                     &
    xc_f90_lda_fxc,                     &
    xc_f90_lda_kxc,                     &
    xc_f90_lda_end,                     &
    xc_f90_lca_init,                    &
    xc_f90_lca_end,                     &
    xc_f90_lca,                         &
    xc_f90_gga_init,                    &
    xc_f90_gga,                         &
    xc_f90_gga_end,                     &
    xc_f90_gga_lb_set_params,           &
    xc_f90_gga_lb_modified,             &
    xc_f90_mgga_init,                   &
    xc_f90_mgga,                        &
    xc_f90_mgga_end

  ! Families of xc functionals
  integer, public, parameter ::     &
    XC_FAMILY_UNKNOWN       =  -1,  &
    XC_FAMILY_LDA           =   1,  &
    XC_FAMILY_GGA           =   2,  &
    XC_FAMILY_MGGA          =   4,  &
    XC_FAMILY_LCA           =   8,  &
    XC_FAMILY_OEP           =  16

  integer, public, parameter ::     &
    XC_UNPOLARIZED          =   1,  &  ! Spin unpolarized
    XC_POLARIZED            =   2      ! Spin polarized

  integer, public, parameter ::     &
    XC_NON_RELATIVISTIC     =   0,  &  ! Functional includes or not realtivistic
    XC_RELATIVISTIC         =   1      ! corrections. Only available in some functionals.

  ! Kinds
  integer, public, parameter ::     &
    XC_EXCHANGE             =   0,  &
    XC_CORRELATION          =   1,  &
    XC_EXCHANGE_CORRELATION =   2


  !----------------------------------------------------------------
  interface
    integer function xc_f90_info_number(info)
      use xc_types_m
      type(xc_info_t), intent(in) :: info
    end function xc_f90_info_number
    
    integer function xc_f90_info_kind(info)
      use xc_types_m
      type(xc_info_t), intent(in) :: info
    end function xc_f90_info_kind

    subroutine xc_f90_info_name(info, s)
      use xc_types_m
      type(xc_info_t),  intent(in)  :: info
      character(len=*), intent(out) :: s
    end subroutine xc_f90_info_name

    integer function xc_f90_info_family(info)
      use xc_types_m
      type(xc_info_t), intent(in)  :: info
    end function xc_f90_info_family

    subroutine xc_f90_info_refs(info, n, s)
      use xc_types_m
      type(xc_info_t),  intent(in)    :: info
      integer,          intent(inout) :: n
      character(len=*), intent(out)   :: s
    end subroutine xc_f90_info_refs
  end interface


  !----------------------------------------------------------------
  interface
    integer function xc_f90_family_from_id(id)
      use xc_types_m
      integer, intent(in) :: id
    end function xc_f90_family_from_id
  end interface


  ! LDAs
  ! We will use the same public interface (xc_lda_init) for the four C procedures
  !----------------------------------------------------------------
  interface xc_f90_lda_init
    subroutine xc_f90_lda_init_(p, info, functional, nspin)
      use xc_types_m
      type(xc_func_t), intent(out) :: p
      type(xc_info_t), intent(out) :: info
      integer,         intent(in)  :: functional
      integer,         intent(in)  :: nspin
    end subroutine xc_f90_lda_init_

    subroutine xc_f90_lda_x_init(p, info, functional, nspin, dim, irel)
      use xc_types_m
      type(xc_func_t), intent(out) :: p
      type(xc_info_t), intent(out) :: info
      integer,         intent(in)  :: functional
      integer,         intent(in)  :: nspin  ! XC_UNPOLARIZED or XC_POLARIZED
      integer,         intent(in)  :: dim    ! 2 or 3 dimensions
      integer,         intent(in)  :: irel   ! XC_NON_RELATIVISTIC or XC_RELATIVISTIC
    end subroutine xc_f90_lda_x_init

    subroutine xc_f90_lda_c_xalpha_init(p, info, functional, nspin, dim, alpha)
      use xc_types_m
      type(xc_func_t), intent(out) :: p
      type(xc_info_t), intent(out) :: info
      integer,         intent(in)  :: functional
      integer,         intent(in)  :: nspin  ! XC_UNPOLARIZED or XC_POLARIZED
      integer,         intent(in)  :: dim    ! 2 or 3 dimensions
      real(8),         intent(in)  :: alpha  ! Ec = alpha Ex
    end subroutine xc_f90_lda_c_xalpha_init

    subroutine xc_f90_lda_c_xalpha_init_sp(p, info, functional, nspin, dim, alpha)
      use xc_types_m
      type(xc_func_t), intent(out) :: p
      type(xc_info_t), intent(out) :: info
      integer,         intent(in)  :: functional
      integer,         intent(in)  :: nspin  ! XC_UNPOLARIZED or XC_POLARIZED
      integer,         intent(in)  :: dim    ! 2 or 3 dimensions
      real(4),         intent(in)  :: alpha  ! Ec = alpha Ex
    end subroutine xc_f90_lda_c_xalpha_init_sp
  end interface


  !----------------------------------------------------------------
  interface
    subroutine xc_f90_lda_end(p)
      use xc_types_m
      type(xc_func_t), intent(inout) :: p
    end subroutine xc_f90_lda_end
  end interface


  !----------------------------------------------------------------
  interface xc_f90_lda
    subroutine xc_f90_lda_dp(p, rho, e, v, fxc, kxc)
      use xc_types_m
      type(xc_func_t), intent(in)  :: p
      real(8),         intent(in)  :: rho   ! rho(nspin) the density
      real(8),         intent(out) :: e     ! the energy per unit particle
      real(8),         intent(out) :: v     ! v(nspin) the potential
      real(8),         intent(out) :: fxc   ! v(nspin,nspin) the xc kernel
      real(8),         intent(out) :: kxc   ! v(nspin,nspin,nspin) the derivative of xc kernel
    end subroutine xc_f90_lda_dp

    subroutine xc_f90_lda_sp(p, rho, e, v, fxc, kxc)
      use xc_types_m
      type(xc_func_t), intent(in)  :: p
      real(4),         intent(in)  :: rho   ! rho(nspin) the density
      real(4),         intent(out) :: e     ! the energy per unit particle
      real(4),         intent(out) :: v     ! v(nspin) the potential
      real(4),         intent(out) :: fxc   ! v(nspin,nspin) the xc kernel
      real(4),         intent(out) :: kxc   ! v(nspin,nspin,nspin) the derivative of xc kernel
    end subroutine xc_f90_lda_sp
  end interface


  !----------------------------------------------------------------
  interface xc_f90_lda_vxc    
    subroutine xc_f90_lda_vxc_dp(p, rho, e, v)
      use xc_types_m
      type(xc_func_t), intent(in)  :: p
      real(8),         intent(in)  :: rho   ! rho(nspin) the density
      real(8),         intent(out) :: e     ! the energy per unit particle
      real(8),         intent(out) :: v     ! v(nspin) the potential
    end subroutine xc_f90_lda_vxc_dp

    subroutine xc_f90_lda_vxc_sp(p, rho, e, v)
      use xc_types_m
      type(xc_func_t), intent(in)  :: p
      real(4),         intent(in)  :: rho   ! rho(nspin) the density
      real(4),         intent(out) :: e     ! the energy per unit particle
      real(4),         intent(out) :: v     ! v(nspin) the potential
    end subroutine xc_f90_lda_vxc_sp
  end interface


  !----------------------------------------------------------------
  interface xc_f90_lda_fxc
    subroutine xc_f90_lda_fxc_dp(p, rho, fxc)
      use xc_types_m
      type(xc_func_t), intent(in)  :: p
      real(8),       intent(in)  :: rho   ! rho(nspin) the density
      real(8),       intent(out) :: fxc   ! v(nspin,nspin) the xc kernel
    end subroutine xc_f90_lda_fxc_dp

    subroutine xc_f90_lda_fxc_sp(p, rho, fxc)
      use xc_types_m
      type(xc_func_t), intent(in)  :: p
      real(4),         intent(in)  :: rho   ! rho(nspin) the density
      real(4),         intent(out) :: fxc   ! v(nspin,nspin) the xc kernel
    end subroutine xc_f90_lda_fxc_sp
  end interface
  

  !----------------------------------------------------------------
  interface xc_f90_lda_kxc
    subroutine xc_f90_lda_kxc_dp(p, rho, kxc)
      use xc_types_m
      type(xc_func_t), intent(in)  :: p
      real(8),         intent(in)  :: rho   ! rho(nspin) the density
      real(8),         intent(out) :: kxc
    end subroutine xc_f90_lda_kxc_dp

    subroutine xc_f90_lda_kxc_sp(p, rho, kxc)
      use xc_types_m
      type(xc_func_t), intent(in)  :: p
      real(4),         intent(in)  :: rho   ! rho(nspin) the density
      real(4),         intent(out) :: kxc
    end subroutine xc_f90_lda_kxc_sp

  end interface
  

  ! GGAs
  ! We will use the same public procedure for the two C procedures.
  !----------------------------------------------------------------
  interface xc_f90_gga_init
    subroutine xc_f90_gga_init_(p, info, functional, nspin)
      use xc_types_m
      type(xc_func_t), intent(out) :: p
      type(xc_info_t), intent(out) :: info
      integer,         intent(in)  :: functional
      integer,         intent(in)  :: nspin
    end subroutine xc_f90_gga_init_
  end interface


  !----------------------------------------------------------------
  interface
    subroutine xc_f90_gga_end(p)
      use xc_types_m
      type(xc_func_t), intent(inout) :: p
    end subroutine xc_f90_gga_end
  end interface


  !----------------------------------------------------------------
  interface xc_f90_gga
    subroutine xc_f90_gga_dp(p, rho, grho, e, dedd, dedgd)
      use xc_types_m
      type(xc_func_t), intent(in)  :: p
      real(8),         intent(in)  :: rho   ! rho(nspin) the density
      real(8),         intent(in)  :: grho  ! grho(3,nspin) the gradient of the density
      real(8),         intent(out) :: e     ! the energy per unit particle
      real(8),         intent(out) :: dedd  ! dedd(nspin) the derivative of the energy
                                            ! in terms of the density
      real(8),         intent(out) :: dedgd ! and in terms of the gradient of the density
    end subroutine xc_f90_gga_dp

    subroutine xc_f90_gga_sp(p, rho, grho, e, dedd, dedgd)
      use xc_types_m
      type(xc_func_t), intent(in)  :: p
      real(4),         intent(in)  :: rho   ! rho(nspin) the density
      real(4),         intent(in)  :: grho  ! grho(3,nspin) the gradient of the density
      real(4),         intent(out) :: e     ! the energy per unit particle
      real(4),         intent(out) :: dedd  ! dedd(nspin) the derivative of the energy
                                          ! in terms of the density
      real(4),         intent(out) :: dedgd ! and in terms of the gradient of the density
    end subroutine xc_f90_gga_sp
  end interface


  !----------------------------------------------------------------
  interface xc_f90_gga_lb_set_params
    subroutine xc_f90_gga_lb_set_params_dp(p, modified, threshold, ip, qtot)
      use xc_types_m
      type(xc_func_t), intent(in)  :: p
      integer,         intent(in)  :: modified   ! should we use the modified version
      real(8),         intent(in)  :: threshold  ! if so, the threshold to use the asymtotic version
      real(8),         intent(in)  :: ip         ! ionization potential
      real(8),         intent(in)  :: qtot       ! total charge
    end subroutine xc_f90_gga_lb_set_params_dp

    subroutine xc_f90_gga_lb_set_params_sp(p, modified, threshold, ip, qtot)
      use xc_types_m
      type(xc_func_t), intent(in)  :: p
      integer,         intent(in)  :: modified   ! should we use the modified version
      real(4),         intent(in)  :: threshold  ! if so, the threshold to use the asymtotic version
      real(4),         intent(in)  :: ip         ! ionization potential
      real(4),         intent(in)  :: qtot       ! total charge
    end subroutine xc_f90_gga_lb_set_params_sp
  end interface


  !----------------------------------------------------------------
  interface xc_f90_gga_lb_modified
    subroutine xc_f90_gga_lb_modified_dp(p, rho,  grho, r, dedd)
      use xc_types_m
      type(xc_func_t), intent(in)  :: p
      real(8),         intent(in)  :: rho   ! rho(nspin) the density
      real(8),         intent(in)  :: grho  ! grho(3,nspin) the gradient of the density
      real(8),         intent(in)  :: r     ! distance from center of finite system
      real(8),         intent(out) :: dedd
    end subroutine xc_f90_gga_lb_modified_dp

    subroutine xc_f90_gga_lb_modified_sp(p, rho,  grho, r, dedd)
      use xc_types_m
      type(xc_func_t), intent(in)  :: p
      real(4),         intent(in)  :: rho   ! rho(nspin) the density
      real(4),         intent(in)  :: grho  ! grho(3,nspin) the gradient of the density
      real(4),         intent(in)  :: r     ! distance from center of finite system
      real(4),         intent(out) :: dedd
    end subroutine xc_f90_gga_lb_modified_sp
  end interface

  ! the meta-GGAs
  !----------------------------------------------------------------
  interface
    subroutine xc_f90_mgga_init(p, info, functional, nspin)
      use xc_types_m
      type(xc_func_t), intent(out) :: p
      type(xc_info_t), intent(out) :: info
      integer,       intent(in)  :: functional
      integer,       intent(in)  :: nspin
    end subroutine xc_f90_mgga_init

    subroutine xc_f90_mgga_end(p)
      use xc_types_m
      type(xc_func_t), intent(inout) :: p
    end subroutine xc_f90_mgga_end
  end interface


  !----------------------------------------------------------------
  interface xc_f90_mgga
    subroutine xc_f90_mgga_dp(p, rho, grho, tau, e, dedd, dedgd, dedtau)
      use xc_types_m
      type(xc_func_t), intent(in)  :: p
      real(8),       intent(in)  :: rho   ! rho(nspin) the density
      real(8),       intent(in)  :: grho  ! grho(3,nspin) the gradient of the density
      real(8),       intent(in)  :: tau   ! tau(nspin) the kinetic energy density
      real(8),       intent(out) :: e     ! the energy per unit particle
      real(8),       intent(out) :: dedd  ! dedd(nspin) the derivative of the energy
      ! in terms of the density
      real(8),       intent(out) :: dedgd ! in terms of the gradient of the density
      real(8),       intent(out) :: dedtau! and in terms of tau
    end subroutine xc_f90_mgga_dp

    subroutine xc_f90_mgga_sp(p, rho, grho, tau, e, dedd, dedgd, dedtau)
      use xc_types_m
      type(xc_func_t), intent(in)  :: p
      real(4),       intent(in)  :: rho   ! rho(nspin) the density
      real(4),       intent(in)  :: grho  ! grho(3,nspin) the gradient of the density
      real(4),       intent(in)  :: tau   ! tau(nspin) the kinetic energy density
      real(4),       intent(out) :: e     ! the energy per unit particle
      real(4),       intent(out) :: dedd  ! dedd(nspin) the derivative of the energy
                                          ! in terms of the density
      real(4),       intent(out) :: dedgd ! in terms of the gradient of the density
      real(4),       intent(out) :: dedtau! and in terms of tau
    end subroutine xc_f90_mgga_sp
  end interface


  ! the LCAs
  !----------------------------------------------------------------
  interface
    subroutine xc_f90_lca_init(p, info, functional, nspin)
      use xc_types_m
      type(xc_func_t), intent(out) :: p
      type(xc_info_t), intent(out) :: info
      integer,       intent(in)  :: functional
      integer,       intent(in)  :: nspin
    end subroutine xc_f90_lca_init

    subroutine xc_f90_lca_end(p)
      use xc_types_m
      type(xc_func_t), intent(inout) :: p
    end subroutine xc_f90_lca_end
  end interface


  !----------------------------------------------------------------
  interface xc_f90_lca
    subroutine xc_f90_lca_dp(p, rho, v, e, dedd, dedv)
      use xc_types_m
      type(xc_func_t), intent(in)  :: p
      real(8),       intent(in)  :: rho   ! rho(nspin) the density
      real(8),       intent(in)  :: v     ! v(3,nspin) the vorticity
      real(8),       intent(out) :: e     ! the energy per unit particle
      real(8),       intent(out) :: dedd  ! dedd(nspin) the derivative of the energy
                                          ! in terms of the density
      real(8),       intent(out) :: dedv  ! and in terms of the vorticity
    end subroutine xc_f90_lca_dp

    subroutine xc_f90_lca_sp(p, rho, v, e, dedd, dedv)
      use xc_types_m
      type(xc_func_t), intent(in)  :: p
      real(4),       intent(in)  :: rho   ! rho(nspin) the density
      real(4),       intent(in)  :: v     ! v(3,nspin) the vorticity
      real(4),       intent(out) :: e     ! the energy per unit particle
      real(4),       intent(out) :: dedd  ! dedd(nspin) the derivative of the energy
                                          ! in terms of the density
      real(4),       intent(out) :: dedv  ! and in terms of the vorticity
    end subroutine xc_f90_lca_sp
  end interface

end module libxc_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
