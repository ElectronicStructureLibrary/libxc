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

module xc_types

  type :: xc_func
     private
     integer, pointer :: func
  end type xc_func

  type :: xc_info
     private
     integer, pointer :: info
  end type xc_info

end module xc_types

module libxc
  
  use xc_types

  implicit none

  private
  public ::                         &
    xc_func,                            &
    xc_info,                            &
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
    xc_f90_gga_lb,                      &
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

  ! the LDAs
  integer, public, parameter ::     &
    XC_LDA_X                =   1,  &  ! Exchange
    XC_LDA_C_WIGNER         =   2,  &  ! Wigner parametrization
    XC_LDA_C_RPA            =   3,  &  ! Random Phase Approximation
    XC_LDA_C_HL             =   4,  &  ! Hedin & Lundqvist
    XC_LDA_C_GL             =   5,  &  ! Gunnarson & Lundqvist
    XC_LDA_C_XALPHA         =   6,  &  ! Slaters Xalpha
    XC_LDA_C_VWN            =   7,  &  ! Vosko, Wilk, & Nussair
    XC_LDA_C_VWN_RPA        =   8,  &  ! Vosko, Wilk, & Nussair (RPA)
    XC_LDA_C_PZ             =   9,  &  ! Perdew & Zunger
    XC_LDA_C_PZ_MOD         =  10,  &  ! Perdew & Zunger (Modified)
    XC_LDA_C_OB_PZ          =  11,  &  ! Ortiz & Ballone (PZ)
    XC_LDA_C_PW             =  12,  &  ! Perdew & Wang
    XC_LDA_C_PW_MOD         =  13,  &  ! Perdew & Wang (Modified)
    XC_LDA_C_OB_PW          =  14,  &  ! Ortiz & Ballone (PW)
    XC_LDA_C_AMGB           =  15      ! Attacalite et al

  ! the GGAs
  integer, public, parameter ::     &
    XC_GGA_X_PBE            = 101,  &  ! Perdew, Burke & Ernzerhof exchange
    XC_GGA_X_PBE_R          = 102,  &  ! Perdew, Burke & Ernzerhof exchange (revised)
    XC_GGA_X_PBE_SOL        = 116,  &  ! Perdew, Burke & Ernzerhof exchange (solids)
    XC_GGA_X_RPBE           = 117,  &  ! Hammer, Hansen & Norskov (PBE-like)
    XC_GGA_X_B86            = 103,  &  ! Becke 86 Xalpha,beta,gamma
    XC_GGA_X_B86_R          = 104,  &  ! Becke 86 Xalpha,beta,gamma reoptimized
    XC_GGA_X_B86_MGC        = 105,  &  ! Becke 88 Xalfa,beta,gamma (with mod. grad. correction)
    XC_GGA_X_B88            = 106,  &  ! Becke 88
    XC_GGA_X_G96            = 107,  &  ! Gill 96
    XC_GGA_X_PW86           = 108,  &  ! Perdew & Wang 86
    XC_GGA_X_PW91           = 109,  &  ! Perdew & Wang 91
    XC_GGA_X_OPTX           = 110,  &  ! Handy & Cohen OPTX 01
    XC_GGA_X_DK87_R1        = 111,  &  ! dePristo & Kress 87 version R1
    XC_GGA_X_DK87_R2        = 112,  &  ! dePristo & Kress 87 version R1
    XC_GGA_X_LG93           = 113,  &  ! Lacks & Gordon 93
    XC_GGA_X_FT97_A         = 114,  &  ! Filatov & Thiel 97 (version A)
    XC_GGA_X_FT97_B         = 115,  &  ! Filatov & Thiel 97 (version B)
    XC_GGA_C_PBE            = 130,  &  ! Perdew, Burke & Ernzerhof correlation
    XC_GGA_C_PBE_SOL        = 133,  &  ! Perdew, Burke & Ernzerhof correlation
    XC_GGA_C_LYP            = 131,  &  ! Lee, Yang & Parr
    XC_GGA_C_P86            = 132,  &  ! Perdew 86
    XC_GGA_XC_LB            = 160      ! van Leeuwen & Baerends


  ! the meta-GGAs
  integer, public, parameter ::     &
    XC_MGGA_X_TPSS          = 201,  &  ! Perdew, Tao, Staroverov & Scuseria exchange
    XC_MGGA_C_TPSS          = 202      ! Perdew, Tao, Staroverov & Scuseria correlation

  ! the LCAs
  integer, public, parameter ::     &
    XC_LCA_OMC              = 301,  &  ! Orestes, Marcasso & Capelle
    XC_LCA_LCH              = 302      ! Lee, Colwell & Handy

  ! info
  interface
    integer function xc_f90_info_number(info)
      use xc_types
      type(xc_info), intent(in) :: info
    end function xc_f90_info_number
    
    integer function xc_f90_info_kind(info)
      use xc_types
      type(xc_info), intent(in) :: info
    end function xc_f90_info_kind

    subroutine xc_f90_info_name(info, s)
      use xc_types
      type(xc_info),         intent(in)  :: info
      character(len=*),      intent(out) :: s
    end subroutine xc_f90_info_name

    integer function xc_f90_info_family(info)
      use xc_types
      type(xc_info), intent(in)  :: info
    end function xc_f90_info_family

    subroutine xc_f90_info_refs(info, n, s)
      use xc_types
      type(xc_info),         intent(in)    :: info
      integer,               intent(inout) :: n
      character(len=*),      intent(out)   :: s
    end subroutine xc_f90_info_refs
  end interface


  ! functionals
  interface
    integer function xc_f90_family_from_id(id)
      use xc_types
      integer, intent(in) :: id
    end function xc_f90_family_from_id
  end interface

  ! LDAs
  ! We will use the same public interface (xc_lda_init) for the four C procedures
  interface xc_f90_lda_init
    subroutine xc_f90_lda_init_(p, info, functional, nspin)
      use xc_types
      type(xc_func),         intent(out) :: p
      type(xc_info),         intent(out) :: info
      integer,               intent(in)  :: functional
      integer,               intent(in)  :: nspin
    end subroutine xc_f90_lda_init_

    subroutine xc_f90_lda_x_init(p, info, functional, nspin, dim, irel)
      use xc_types
      type(xc_func), intent(out) :: p
      type(xc_info), intent(out) :: info
      integer,       intent(in)  :: functional
      integer,       intent(in)  :: nspin  ! XC_UNPOLARIZED or XC_POLARIZED
      integer,       intent(in)  :: dim    ! 2 or 3 dimensions
      integer,       intent(in)  :: irel   ! XC_NON_RELATIVISTIC or XC_RELATIVISTIC
    end subroutine xc_f90_lda_x_init

    subroutine xc_f90_lda_c_xalpha_init(p, info, functional, nspin, dim, alpha)
      use xc_types
      type(xc_func), intent(out) :: p
      type(xc_info), intent(out) :: info
      integer,       intent(in)  :: functional
      integer,       intent(in)  :: nspin  ! XC_UNPOLARIZED or XC_POLARIZED
      integer,       intent(in)  :: dim    ! 2 or 3 dimensions
      real(8),       intent(in)  :: alpha  ! Ec = alpha Ex
    end subroutine xc_f90_lda_c_xalpha_init

    subroutine xc_f90_lda_c_xalpha_init_sp(p, info, functional, nspin, dim, alpha)
      use xc_types
      type(xc_func), intent(out) :: p
      type(xc_info), intent(out) :: info
      integer,       intent(in)  :: functional
      integer,       intent(in)  :: nspin  ! XC_UNPOLARIZED or XC_POLARIZED
      integer,       intent(in)  :: dim    ! 2 or 3 dimensions
      real(4),       intent(in)  :: alpha  ! Ec = alpha Ex
    end subroutine xc_f90_lda_c_xalpha_init_sp
  end interface

  interface
    subroutine xc_f90_lda_end(p)
      use xc_types
      type(xc_func), intent(inout) :: p
    end subroutine xc_f90_lda_end
  end interface

  interface xc_f90_lda

    subroutine xc_f90_lda_dp(p, rho, e, v, fxc, kxc)
      use xc_types
      type(xc_func), intent(in)  :: p
      real(8),       intent(in)  :: rho   ! rho(nspin) the density
      real(8),       intent(out) :: e     ! the energy per unit particle
      real(8),       intent(out) :: v     ! v(nspin) the potential
      real(8),       intent(out) :: fxc   ! v(nspin,nspin) the xc kernel
      real(8),       intent(out) :: kxc   ! v(nspin,nspin,nspin) the derivative of xc kernel
    end subroutine xc_f90_lda_dp

    subroutine xc_f90_lda_sp(p, rho, e, v, fxc, kxc)
      use xc_types
      type(xc_func), intent(in)  :: p
      real(4),       intent(in)  :: rho   ! rho(nspin) the density
      real(4),       intent(out) :: e     ! the energy per unit particle
      real(4),       intent(out) :: v     ! v(nspin) the potential
      real(4),       intent(out) :: fxc   ! v(nspin,nspin) the xc kernel
      real(4),       intent(out) :: kxc   ! v(nspin,nspin,nspin) the derivative of xc kernel
    end subroutine xc_f90_lda_sp

  end interface

  interface xc_f90_lda_vxc
    
    subroutine xc_f90_lda_vxc_dp(p, rho, e, v)
      use xc_types
      type(xc_func), intent(in)  :: p
      real(8),       intent(in)  :: rho   ! rho(nspin) the density
      real(8),       intent(out) :: e     ! the energy per unit particle
      real(8),       intent(out) :: v     ! v(nspin) the potential
    end subroutine xc_f90_lda_vxc_dp

    subroutine xc_f90_lda_vxc_sp(p, rho, e, v)
      use xc_types
      type(xc_func), intent(in)  :: p
      real(4),       intent(in)  :: rho   ! rho(nspin) the density
      real(4),       intent(out) :: e     ! the energy per unit particle
      real(4),       intent(out) :: v     ! v(nspin) the potential
    end subroutine xc_f90_lda_vxc_sp

  end interface

  interface xc_f90_lda_fxc

    subroutine xc_f90_lda_fxc_dp(p, rho, fxc)
      use xc_types
      type(xc_func), intent(in)  :: p
      real(8),       intent(in)  :: rho   ! rho(nspin) the density
      real(8),       intent(out) :: fxc   ! v(nspin,nspin) the xc kernel
    end subroutine xc_f90_lda_fxc_dp

    subroutine xc_f90_lda_fxc_sp(p, rho, fxc)
      use xc_types
      type(xc_func), intent(in)  :: p
      real(4),       intent(in)  :: rho   ! rho(nspin) the density
      real(4),       intent(out) :: fxc   ! v(nspin,nspin) the xc kernel
    end subroutine xc_f90_lda_fxc_sp

  end interface
  
  interface xc_f90_lda_kxc

    subroutine xc_f90_lda_kxc_dp(p, rho, kxc)
      use xc_types
      type(xc_func), intent(in)  :: p
      real(8),       intent(in)  :: rho   ! rho(nspin) the density
      real(8),       intent(out) :: kxc
    end subroutine xc_f90_lda_kxc_dp

    subroutine xc_f90_lda_kxc_sp(p, rho, kxc)
      use xc_types
      type(xc_func), intent(in)  :: p
      real(4),       intent(in)  :: rho   ! rho(nspin) the density
      real(4),       intent(out) :: kxc
    end subroutine xc_f90_lda_kxc_sp

  end interface
  
  ! GGAs
  ! We will use the same public procedure for the two C procedures.
  interface xc_f90_gga_init
    subroutine xc_f90_gga_init_(p, info, functional, nspin)
      use xc_types
      type(xc_func), intent(out) :: p
      type(xc_info), intent(out) :: info
      integer,   intent(in)  :: functional
      integer,   intent(in)  :: nspin
    end subroutine xc_f90_gga_init_

    subroutine xc_f90_gga_lb_init(p, info, functional, nspin, modified, threshold)
      use xc_types
      type(xc_func), intent(out) :: p
      type(xc_info), intent(out) :: info
      integer,       intent(in)  :: functional
      integer,       intent(in)  :: nspin
      integer,       intent(in)  :: modified
      real(8),       intent(in)  :: threshold
    end subroutine xc_f90_gga_lb_init

    subroutine xc_f90_gga_lb_init_sp(p, info, functional, nspin, modified, threshold)
      use xc_types
      type(xc_func), intent(out) :: p
      type(xc_info), intent(out) :: info
      integer,       intent(in)  :: functional
      integer,       intent(in)  :: nspin
      integer,       intent(in)  :: modified
      real(4),       intent(in)  :: threshold
    end subroutine xc_f90_gga_lb_init_sp
  end interface

  interface
    subroutine xc_f90_gga_end(p)
      use xc_types
      type(xc_func), intent(inout) :: p
    end subroutine xc_f90_gga_end
  end interface

  interface xc_f90_gga

    subroutine xc_f90_gga_dp(p, rho, grho, e, dedd, dedgd)
      use xc_types
      type(xc_func), intent(in)  :: p
      real(8),       intent(in)  :: rho   ! rho(nspin) the density
      real(8),       intent(in)  :: grho  ! grho(3,nspin) the gradient of the density
      real(8),       intent(out) :: e     ! the energy per unit particle
      real(8),       intent(out) :: dedd  ! dedd(nspin) the derivative of the energy
                                          ! in terms of the density
      real(8),       intent(out) :: dedgd ! and in terms of the gradient of the density
    end subroutine xc_f90_gga_dp

    subroutine xc_f90_gga_sp(p, rho, grho, e, dedd, dedgd)
      use xc_types
      type(xc_func), intent(in)  :: p
      real(4),       intent(in)  :: rho   ! rho(nspin) the density
      real(4),       intent(in)  :: grho  ! grho(3,nspin) the gradient of the density
      real(4),       intent(out) :: e     ! the energy per unit particle
      real(4),       intent(out) :: dedd  ! dedd(nspin) the derivative of the energy
                                          ! in terms of the density
      real(4),       intent(out) :: dedgd ! and in terms of the gradient of the density
    end subroutine xc_f90_gga_sp

  end interface

  interface xc_f90_gga_lb

    subroutine xc_f90_gga_lb_dp(p, rho, grho, r, ip, qtot, dedd)
      use xc_types
      type(xc_func), intent(in)  :: p
      real(8),       intent(in)  :: rho   ! rho(nspin) the density
      real(8),       intent(in)  :: grho  ! grho(3,nspin) the gradient of the density
      real(8),       intent(in)  :: r     ! distance from center of finite system
      real(8),       intent(in)  :: ip    ! ionization potential
      real(8),       intent(in)  :: qtot  ! total charge
      real(8),       intent(out) :: dedd
    end subroutine xc_f90_gga_lb_dp

    subroutine xc_f90_gga_lb_sp(p, rho, grho, r, ip, qtot, dedd)
      use xc_types
      type(xc_func), intent(in)  :: p
      real(4),       intent(in)  :: rho   ! rho(nspin) the density
      real(4),       intent(in)  :: grho  ! grho(3,nspin) the gradient of the density
      real(4),       intent(in)  :: r     ! distance from center of finite system
      real(4),       intent(in)  :: ip    ! ionization potential
      real(4),       intent(in)  :: qtot  ! total charge
      real(4),       intent(out) :: dedd
    end subroutine xc_f90_gga_lb_sp

  end interface


  ! the meta-GGAs
  interface
    subroutine xc_f90_mgga_init(p, info, functional, nspin)
      use xc_types
      type(xc_func), intent(out) :: p
      type(xc_info), intent(out) :: info
      integer,       intent(in)  :: functional
      integer,       intent(in)  :: nspin
    end subroutine xc_f90_mgga_init

    subroutine xc_f90_mgga_end(p)
      use xc_types
      type(xc_func), intent(inout) :: p
    end subroutine xc_f90_mgga_end
  end interface

  interface xc_f90_mgga

    subroutine xc_f90_mgga_dp(p, rho, grho, tau, e, dedd, dedgd, dedtau)
      use xc_types
      type(xc_func), intent(in)  :: p
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
      use xc_types
      type(xc_func), intent(in)  :: p
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
  interface
    subroutine xc_f90_lca_init(p, info, functional, nspin)
      use xc_types
      type(xc_func), intent(out) :: p
      type(xc_info), intent(out) :: info
      integer,       intent(in)  :: functional
      integer,       intent(in)  :: nspin
    end subroutine xc_f90_lca_init

    subroutine xc_f90_lca_end(p)
      use xc_types
      type(xc_func), intent(inout) :: p
    end subroutine xc_f90_lca_end
  end interface

  interface xc_f90_lca

    subroutine xc_f90_lca_dp(p, rho, v, e, dedd, dedv)
      use xc_types
      type(xc_func), intent(in)  :: p
      real(8),       intent(in)  :: rho   ! rho(nspin) the density
      real(8),       intent(in)  :: v     ! v(3,nspin) the vorticity
      real(8),       intent(out) :: e     ! the energy per unit particle
      real(8),       intent(out) :: dedd  ! dedd(nspin) the derivative of the energy
                                          ! in terms of the density
      real(8),       intent(out) :: dedv  ! and in terms of the vorticity
    end subroutine xc_f90_lca_dp

    subroutine xc_f90_lca_sp(p, rho, v, e, dedd, dedv)
      use xc_types
      type(xc_func), intent(in)  :: p
      real(4),       intent(in)  :: rho   ! rho(nspin) the density
      real(4),       intent(in)  :: v     ! v(3,nspin) the vorticity
      real(4),       intent(out) :: e     ! the energy per unit particle
      real(4),       intent(out) :: dedd  ! dedd(nspin) the derivative of the energy
                                          ! in terms of the density
      real(4),       intent(out) :: dedv  ! and in terms of the vorticity
    end subroutine xc_f90_lca_sp

  end interface

end module libxc

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
