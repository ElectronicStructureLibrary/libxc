!! Copyright (C) 2003-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
!!
!! This program is free software; you can redistribute it and/or modify
!! it under the terms of the GNU Lesser General Public License as published by
!! the Free Software Foundation; either version 2, or (at your option)
!! any later version.
!!
!! This program is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU Lesser General Public License for more details.
!!
!! You should have received a copy of the GNU Lesser General Public License
!! along with this program; if not, write to the Free Software
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.
!!
!! $Id: libxc.f90 3550 2007-11-19 14:32:49Z marques $

#if SINGLE_PRECISION
#  define XC_F90(x) xc_s_f90_ ## x
#else
#  define XC_F90(x) xc_f90_ ## x
#endif

module XC_F90(types_m)
#if SINGLE_PRECISION
  integer, public, parameter :: xc_f90_kind = selected_real_kind(4)
#else
  integer, public, parameter :: xc_f90_kind = selected_real_kind(14)
#endif

  type XC_F90(func_t)
    private
    integer, pointer :: func
  end type XC_F90(func_t)

  type XC_F90(info_t)
    private
    integer, pointer :: info
  end type XC_F90(info_t)

end module XC_F90(types_m)

module XC_F90(lib_m)
  
  use XC_F90(types_m)
  use libxc_funcs_m

  implicit none

  public
  public ::                              &
    XC_F90(func_t),                      &
    XC_F90(info_t),                      &
    XC_F90(info_number),                 &
    XC_F90(info_kind),                   &
    XC_F90(info_name),                   &
    XC_F90(info_family),                 &
    XC_F90(info_provides),               &
    XC_F90(info_refs),                   &
    XC_F90(family_from_id),              &
    XC_F90(lda_init),                    &
    XC_F90(lda),                         &
    XC_F90(lda_exc),                     &
    XC_F90(lda_exc_vxc),                 &
    XC_F90(lda_vxc),                     &
    XC_F90(lda_fxc),                     &
    XC_F90(lda_kxc),                     &
    XC_F90(lda_end),                     &
    XC_F90(lda_c_1d_csc_set_par),        &
    XC_F90(lda_c_2d_prm_set_par),        &
    XC_F90(gga_init),                    &
    XC_F90(gga),                         &
    XC_F90(gga_exc),                     &
    XC_F90(gga_exc_vxc),                 &
    XC_F90(gga_vxc),                     &
    XC_F90(gga_fxc),                     &
    XC_F90(gga_end),                     &
    XC_F90(gga_lb_set_par),              &
    XC_F90(gga_lb_modified),             &
    XC_F90(hyb_gga_init),                &
    XC_F90(hyb_gga_end),                 &
    XC_F90(hyb_gga),                     &
    XC_F90(hyb_gga_exc),                 &
    XC_F90(hyb_gga_exc_vxc),             &
    XC_F90(hyb_gga_vxc),                 &
    XC_F90(hyb_gga_fxc),                 &
    XC_F90(hyb_gga_exx_coef),            &
    XC_F90(lca_init),                    &
    XC_F90(lca_end),                     &
    XC_F90(lca),                         &
    XC_F90(mgga_init),                   &
    XC_F90(mgga),                        &
    XC_F90(mgga_exc),                    &
    XC_F90(mgga_exc_vxc),                    &
    XC_F90(mgga_vxc),                    &
    XC_F90(mgga_fxc),                    &
    XC_F90(mgga_end)

  ! Families of xc functionals
  integer, public, parameter ::     &
    XC_FAMILY_UNKNOWN       =  -1,  &
    XC_FAMILY_NONE          =   0,  &
    XC_FAMILY_LDA           =   1,  &
    XC_FAMILY_GGA           =   2,  &
    XC_FAMILY_MGGA          =   4,  &
    XC_FAMILY_LCA           =   8,  &
    XC_FAMILY_OEP           =  16,  &
    XC_FAMILY_HYB_GGA       =  32

  integer, public, parameter ::     &
    XC_UNPOLARIZED          =   1,  &  ! Spin unpolarized
    XC_POLARIZED            =   2      ! Spin polarized

  integer, public, parameter ::     &
    XC_NON_RELATIVISTIC     =   0,  &  ! Functional includes or not relativistic
    XC_RELATIVISTIC         =   1      ! corrections. Only available in some functionals.

  ! Kinds
  integer, public, parameter ::     &
    XC_EXCHANGE             =   0,  &
    XC_CORRELATION          =   1,  &
    XC_EXCHANGE_CORRELATION =   2

  integer, public, parameter ::     &
    XC_PROVIDES_EXC         =   1,  &
    XC_PROVIDES_VXC         =   2,  &
    XC_PROVIDES_FXC         =   4,  &
    XC_PROVIDES_KXC         =   8
 

  !----------------------------------------------------------------
  interface
    integer function XC_F90(info_number)(info)
      use XC_F90(types_m)
      type(XC_F90(info_t)), intent(in) :: info
    end function XC_F90(info_number)
    
    integer function XC_F90(info_kind)(info)
      use XC_F90(types_m)
      type(XC_F90(info_t)), intent(in) :: info
    end function XC_F90(info_kind)

    subroutine XC_F90(info_name)(info, s)
      use XC_F90(types_m)
      type(XC_F90(info_t)),  intent(in)  :: info
      character(len=*), intent(out) :: s
    end subroutine XC_F90(info_name)

    integer function XC_F90(info_family)(info)
      use XC_F90(types_m)
      type(XC_F90(info_t)), intent(in)  :: info
    end function XC_F90(info_family)

    integer function XC_F90(info_provides)(info)
      use XC_F90(types_m)
      type(XC_F90(info_t)), intent(in)  :: info
    end function XC_F90(info_provides)

    subroutine XC_F90(info_refs)(info, n, s)
      use XC_F90(types_m)
      type(XC_F90(info_t)),  intent(in)    :: info
      integer,          intent(inout) :: n
      character(len=*), intent(out)   :: s
    end subroutine XC_F90(info_refs)

    integer function XC_F90(family_from_id)(id)
      use XC_F90(types_m)
      integer, intent(in) :: id
    end function XC_F90(family_from_id)
  end interface


  ! LDAs
  ! We will use the same public interface (xc_lda_init) for the 3 C procedures
  !----------------------------------------------------------------
  interface XC_F90(lda_init)
    subroutine XC_F90(lda_init_)(p, info, functional, nspin)
      use XC_F90(types_m)
      type(XC_F90(func_t)), intent(out) :: p
      type(XC_F90(info_t)), intent(out) :: info
      integer,         intent(in)  :: functional
      integer,         intent(in)  :: nspin
    end subroutine XC_F90(lda_init_)

    subroutine XC_F90(lda_x_init)(p, info, functional, nspin, dim, irel)
      use XC_F90(types_m)
      type(XC_F90(func_t)), intent(out) :: p
      type(XC_F90(info_t)), intent(out) :: info
      integer,         intent(in)  :: functional
      integer,         intent(in)  :: nspin  ! XC_UNPOLARIZED or XC_POLARIZED
      integer,         intent(in)  :: dim    ! 2 or 3 dimensions
      integer,         intent(in)  :: irel   ! XC_NON_RELATIVISTIC or XC_RELATIVISTIC
    end subroutine XC_F90(lda_x_init)

    subroutine XC_F90(lda_c_xalpha_init)(p, info, functional, nspin, dim, alpha)
      use XC_F90(types_m)
      type(XC_F90(func_t)),   intent(out) :: p
      type(XC_F90(info_t)),   intent(out) :: info
      integer,           intent(in)  :: functional
      integer,           intent(in)  :: nspin  ! XC_UNPOLARIZED or XC_POLARIZED
      integer,           intent(in)  :: dim    ! 2 or 3 dimensions
      real(xc_f90_kind), intent(in)  :: alpha  ! Ec = alpha Ex
    end subroutine XC_F90(lda_c_xalpha_init)
  end interface


  !----------------------------------------------------------------
  interface
    subroutine XC_F90(lda_end)(p)
      use XC_F90(types_m)
      type(XC_F90(func_t)), intent(inout) :: p
    end subroutine XC_F90(lda_end)
  end interface


  !----------------------------------------------------------------
  interface
    subroutine XC_F90(lda)(p, np, rho, zk, vrho, fxc, kxc)
      use XC_F90(types_m)
      type(XC_F90(func_t)), intent(in)  :: p
      integer,              intent(in)  :: np
      real(xc_f90_kind),    intent(in)  :: rho   ! rho(nspin) the density
      real(xc_f90_kind),    intent(out) :: zk    ! the energy per unit particle
      real(xc_f90_kind),    intent(out) :: vrho  ! v(nspin) the potential
      real(xc_f90_kind),    intent(out) :: fxc   ! v(nspin,nspin) the xc kernel
      real(xc_f90_kind),    intent(out) :: kxc   ! v(nspin,nspin,nspin) the derivative of xc kernel
    end subroutine XC_F90(lda)

    subroutine XC_F90(lda_exc)(p, np, rho, zk)
      use XC_F90(types_m)
      type(XC_F90(func_t)), intent(in)  :: p
      integer,              intent(in)  :: np
      real(xc_f90_kind),    intent(in)  :: rho   ! rho(nspin) the density
      real(xc_f90_kind),    intent(out) :: zk    ! the energy per unit particle
    end subroutine XC_F90(lda_exc)

    subroutine XC_F90(lda_exc_vxc)(p, np, rho, e, v)
      use XC_F90(types_m)
      type(XC_F90(func_t)), intent(in)  :: p
      integer,              intent(in)  :: np
      real(xc_f90_kind),    intent(in)  :: rho   ! rho(nspin) the density
      real(xc_f90_kind),    intent(out) :: e     ! the energy per unit particle
      real(xc_f90_kind),    intent(out) :: v     ! v(nspin) the potential
    end subroutine XC_F90(lda_exc_vxc)

    subroutine XC_F90(lda_vxc)(p, np, rho, v)
      use XC_F90(types_m)
      type(XC_F90(func_t)), intent(in)  :: p
      integer,              intent(in)  :: np
      real(xc_f90_kind),    intent(in)  :: rho   ! rho(nspin) the density
      real(xc_f90_kind),    intent(out) :: v     ! v(nspin) the potential
    end subroutine XC_F90(lda_vxc)

    subroutine XC_F90(lda_fxc)(p, np, rho, fxc)
      use XC_F90(types_m)
      type(XC_F90(func_t)), intent(in)  :: p
      integer,              intent(in)  :: np
      real(xc_f90_kind),    intent(in)  :: rho   ! rho(nspin) the density
      real(xc_f90_kind),    intent(out) :: fxc   ! v(nspin,nspin) the xc kernel
    end subroutine XC_F90(lda_fxc)

    subroutine XC_F90(lda_kxc)(p, np, rho, kxc)
      use XC_F90(types_m)
      type(XC_F90(func_t)), intent(in)  :: p
      integer,              intent(in)  :: np
      real(xc_f90_kind),    intent(in)  :: rho   ! rho(nspin) the density
      real(xc_f90_kind),    intent(out) :: kxc
    end subroutine XC_F90(lda_kxc)
  end interface
  

  interface
    subroutine XC_F90(lda_c_1d_csc_set_par)(p, bb)
      use XC_F90(types_m)
      type(XC_F90(func_t)), intent(inout)  :: p
      real(xc_f90_kind),    intent(in)     :: bb
    end subroutine XC_F90(lda_c_1d_csc_set_par)
  end interface

  interface
    subroutine XC_F90(lda_c_2d_prm_set_par)(p, N)
      use XC_F90(types_m)
      type(XC_F90(func_t)), intent(inout)  :: p
      real(xc_f90_kind),    intent(in)     :: N
    end subroutine XC_F90(lda_c_2d_prm_set_par)
  end interface

  ! GGAs
  !----------------------------------------------------------------
  interface
    subroutine XC_F90(gga_init)(p, info, functional, nspin)
      use XC_F90(types_m)
      type(XC_F90(func_t)), intent(out) :: p
      type(XC_F90(info_t)), intent(out) :: info
      integer,              intent(in)  :: functional
      integer,              intent(in)  :: nspin
    end subroutine XC_F90(gga_init)

    subroutine XC_F90(gga_end)(p)
      use XC_F90(types_m)
      type(XC_F90(func_t)), intent(inout) :: p
    end subroutine XC_F90(gga_end)
  end interface


  !----------------------------------------------------------------
  interface
    subroutine XC_F90(gga)(p, rho, sigma, zk, vrho, vsigma, v2rho2, v2rhosigma, v2sigma2)
      use XC_F90(types_m)
      type(XC_F90(func_t)), intent(in)  :: p
      real(xc_f90_kind),    intent(in)  :: rho
      real(xc_f90_kind),    intent(in)  :: sigma
      real(xc_f90_kind),    intent(out) :: zk
      real(xc_f90_kind),    intent(out) :: vrho
      real(xc_f90_kind),    intent(out) :: vsigma
      real(xc_f90_kind),    intent(out) :: v2rho2
      real(xc_f90_kind),    intent(out) :: v2rhosigma
      real(xc_f90_kind),    intent(out) :: v2sigma2
    end subroutine XC_F90(gga)

    subroutine XC_F90(gga_exc)(p, rho, sigma, zk)
      use XC_F90(types_m)
      type(XC_F90(func_t)), intent(in)  :: p
      real(xc_f90_kind),    intent(in)  :: rho
      real(xc_f90_kind),    intent(in)  :: sigma
      real(xc_f90_kind),    intent(out) :: zk
    end subroutine XC_F90(gga_exc)

    subroutine XC_F90(gga_exc_vxc)(p, rho, sigma, zk, vrho, vsigma)
      use XC_F90(types_m)
      type(XC_F90(func_t)), intent(in)  :: p
      real(xc_f90_kind),    intent(in)  :: rho
      real(xc_f90_kind),    intent(in)  :: sigma
      real(xc_f90_kind),    intent(out) :: zk
      real(xc_f90_kind),    intent(out) :: vrho
      real(xc_f90_kind),    intent(out) :: vsigma
    end subroutine XC_F90(gga_exc_vxc)

    subroutine XC_F90(gga_vxc)(p, rho, sigma, vrho, vsigma)
      use XC_F90(types_m)
      type(XC_F90(func_t)), intent(in)  :: p
      real(xc_f90_kind),    intent(in)  :: rho
      real(xc_f90_kind),    intent(in)  :: sigma
      real(xc_f90_kind),    intent(out) :: vrho
      real(xc_f90_kind),    intent(out) :: vsigma
    end subroutine XC_F90(gga_vxc)

    subroutine XC_F90(gga_fxc)(p, rho, sigma, v2rho2, v2rhosigma, v2sigma2)
      use XC_F90(types_m)
      type(XC_F90(func_t)), intent(in)  :: p
      real(xc_f90_kind),    intent(in)  :: rho
      real(xc_f90_kind),    intent(in)  :: sigma
      real(xc_f90_kind),    intent(out) :: v2rho2
      real(xc_f90_kind),    intent(out) :: v2rhosigma
      real(xc_f90_kind),    intent(out) :: v2sigma2
    end subroutine XC_F90(gga_fxc)
  end interface

  !----------------------------------------------------------------
  interface
    subroutine XC_F90(gga_lb_set_par)(p, modified, threshold, ip, qtot)
      use XC_F90(types_m)
      type(XC_F90(func_t)),   intent(in)  :: p
      integer,           intent(in)  :: modified   ! should we use the modified version
      real(xc_f90_kind), intent(in)  :: threshold  ! if so, the threshold to use the asymtotic version
      real(xc_f90_kind), intent(in)  :: ip         ! ionization potential
      real(xc_f90_kind), intent(in)  :: qtot       ! total charge
    end subroutine XC_F90(gga_lb_set_par)
  end interface


  !----------------------------------------------------------------
  interface
    subroutine XC_F90(gga_lb_modified)(p, rho,  grho, r, dedd)
      use XC_F90(types_m)
      type(XC_F90(func_t)),   intent(in)  :: p
      real(xc_f90_kind), intent(in)  :: rho   ! rho(nspin) the density
      real(xc_f90_kind), intent(in)  :: grho  ! grho(3,nspin) the gradient of the density
      real(xc_f90_kind), intent(in)  :: r     ! distance from center of finite system
      real(xc_f90_kind), intent(out) :: dedd
    end subroutine XC_F90(gga_lb_modified)
  end interface


  ! Hybrids GGAs
  !----------------------------------------------------------------
  interface
    subroutine XC_F90(hyb_gga_init)(p, info, functional, nspin)
      use XC_F90(types_m)
      type(XC_F90(func_t)), intent(out) :: p
      type(XC_F90(info_t)), intent(out) :: info
      integer,         intent(in)  :: functional
      integer,         intent(in)  :: nspin
    end subroutine XC_F90(hyb_gga_init)

    subroutine XC_F90(hyb_gga_end)(p)
      use XC_F90(types_m)
      type(XC_F90(func_t)), intent(inout) :: p
    end subroutine XC_F90(hyb_gga_end)
  end interface


  !----------------------------------------------------------------
  interface
    subroutine XC_F90(hyb_gga)(p, rho, grho, e, dedd, dedgd)
      use XC_F90(types_m)
      type(XC_F90(func_t)),   intent(in)  :: p
      real(xc_f90_kind), intent(in)  :: rho   ! rho(nspin) the density
      real(xc_f90_kind), intent(in)  :: grho  ! grho(3,nspin) the gradient of the density
      real(xc_f90_kind), intent(out) :: e     ! the energy per unit particle
      real(xc_f90_kind), intent(out) :: dedd  ! dedd(nspin) the derivative of the energy
                                            ! in terms of the density
      real(xc_f90_kind), intent(out) :: dedgd ! and in terms of the gradient of the density
    end subroutine XC_F90(hyb_gga)

    subroutine XC_F90(hyb_gga_exc)(p, rho, sigma, zk)
      use XC_F90(types_m)
      type(XC_F90(func_t)), intent(in)  :: p
      real(xc_f90_kind),    intent(in)  :: rho
      real(xc_f90_kind),    intent(in)  :: sigma
      real(xc_f90_kind),    intent(out) :: zk
    end subroutine XC_F90(hyb_gga_exc)

    subroutine XC_F90(hyb_gga_exc_vxc)(p, rho, sigma, zk, vrho, vsigma)
      use XC_F90(types_m)
      type(XC_F90(func_t)), intent(in)  :: p
      real(xc_f90_kind),    intent(in)  :: rho
      real(xc_f90_kind),    intent(in)  :: sigma
      real(xc_f90_kind),    intent(out) :: zk
      real(xc_f90_kind),    intent(out) :: vrho
      real(xc_f90_kind),    intent(out) :: vsigma
    end subroutine XC_F90(hyb_gga_exc_vxc)

    subroutine XC_F90(hyb_gga_vxc)(p, rho, sigma, vrho, vsigma)
      use XC_F90(types_m)
      type(XC_F90(func_t)), intent(in)  :: p
      real(xc_f90_kind),    intent(in)  :: rho
      real(xc_f90_kind),    intent(in)  :: sigma
      real(xc_f90_kind),    intent(out) :: vrho
      real(xc_f90_kind),    intent(out) :: vsigma
    end subroutine XC_F90(hyb_gga_vxc)

    subroutine XC_F90(hyb_gga_fxc)(p, rho, sigma, v2rho2, v2rhosigma, v2sigma2)
      use XC_F90(types_m)
      type(XC_F90(func_t)), intent(in)  :: p
      real(xc_f90_kind),    intent(in)  :: rho
      real(xc_f90_kind),    intent(in)  :: sigma
      real(xc_f90_kind),    intent(out) :: v2rho2
      real(xc_f90_kind),    intent(out) :: v2rhosigma
      real(xc_f90_kind),    intent(out) :: v2sigma2
    end subroutine XC_F90(hyb_gga_fxc)
  end interface


  !----------------------------------------------------------------
  interface
    subroutine XC_F90(hyb_gga_exx_coef)(p, coef)
      use XC_F90(types_m)
      type(XC_F90(func_t)),   intent(in)  :: p
      real(xc_f90_kind), intent(out) :: coef
    end subroutine XC_F90(hyb_gga_exx_coef)
  end interface


  ! the meta-GGAs
  !----------------------------------------------------------------
  interface
    subroutine XC_F90(mgga_init)(p, info, functional, nspin)
      use XC_F90(types_m)
      type(XC_F90(func_t)), intent(out) :: p
      type(XC_F90(info_t)), intent(out) :: info
      integer,              intent(in)  :: functional
      integer,              intent(in)  :: nspin
    end subroutine XC_F90(mgga_init)

    subroutine XC_F90(mgga_end)(p)
      use XC_F90(types_m)
      type(XC_F90(func_t)), intent(inout) :: p
    end subroutine XC_F90(mgga_end)
  end interface


  !----------------------------------------------------------------
  interface
    subroutine XC_F90(mgga)(p, rho, sigma, lrho, tau, zk, vrho, vsigma, vlrho, vtau, &
      v2rho2, v2rhosigma, v2sigma2, v2rhotau, v2tausigma, v2tau2)

      use XC_F90(types_m)
      type(XC_F90(func_t)), intent(in)  :: p
      real(xc_f90_kind),    intent(in)  :: rho   
      real(xc_f90_kind),    intent(in)  :: sigma 
      real(xc_f90_kind),    intent(in)  :: lrho
      real(xc_f90_kind),    intent(in)  :: tau   
      real(xc_f90_kind),    intent(out) :: zk
      real(xc_f90_kind),    intent(out) :: vrho
      real(xc_f90_kind),    intent(out) :: vsigma
      real(xc_f90_kind),    intent(out) :: vlrho
      real(xc_f90_kind),    intent(out) :: vtau
      real(xc_f90_kind),    intent(out) :: v2rho2
      real(xc_f90_kind),    intent(out) :: v2rhosigma
      real(xc_f90_kind),    intent(out) :: v2sigma2
      real(xc_f90_kind),    intent(out) :: v2rhotau
      real(xc_f90_kind),    intent(out) :: v2tausigma
      real(xc_f90_kind),    intent(out) :: v2tau2
    end subroutine XC_F90(mgga)

    subroutine XC_F90(mgga_exc)(p, rho, sigma, lrho, tau, zk)
      use XC_F90(types_m)
      type(XC_F90(func_t)), intent(in)  :: p
      real(xc_f90_kind),    intent(in)  :: rho
      real(xc_f90_kind),    intent(in)  :: sigma
      real(xc_f90_kind),    intent(in)  :: lrho
      real(xc_f90_kind),    intent(in)  :: tau
      real(xc_f90_kind),    intent(out) :: zk
    end subroutine XC_F90(mgga_exc)

    subroutine XC_F90(mgga_exc_vxc)(p, rho, sigma, lrho, tau, zk, vrho, vsigma, vlrho, vtau)
      use XC_F90(types_m)
      type(XC_F90(func_t)), intent(in)  :: p
      real(xc_f90_kind),    intent(in)  :: rho   
      real(xc_f90_kind),    intent(in)  :: sigma 
      real(xc_f90_kind),    intent(in)  :: lrho   
      real(xc_f90_kind),    intent(in)  :: tau
      real(xc_f90_kind),    intent(out) :: zk
      real(xc_f90_kind),    intent(out) :: vrho
      real(xc_f90_kind),    intent(out) :: vsigma
      real(xc_f90_kind),    intent(out) :: vlrho
      real(xc_f90_kind),    intent(out) :: vtau
    end subroutine XC_F90(mgga_exc_vxc)

    subroutine XC_F90(mgga_vxc)(p, rho, sigma, lrho, tau, vrho, vsigma, vlrho, vtau)
      use XC_F90(types_m)
      type(XC_F90(func_t)), intent(in)  :: p
      real(xc_f90_kind),    intent(in)  :: rho   
      real(xc_f90_kind),    intent(in)  :: sigma 
      real(xc_f90_kind),    intent(in)  :: lrho
      real(xc_f90_kind),    intent(in)  :: tau
      real(xc_f90_kind),    intent(out) :: vrho
      real(xc_f90_kind),    intent(out) :: vsigma
      real(xc_f90_kind),    intent(out) :: vlrho
      real(xc_f90_kind),    intent(out) :: vtau
    end subroutine XC_F90(mgga_vxc)

    subroutine XC_F90(mgga_fxc)(p, rho, sigma, lrho, tau, &
      v2rho2, v2rhosigma, v2sigma2, v2rhotau, v2tausigma, v2tau2)

      use XC_F90(types_m)
      type(XC_F90(func_t)), intent(in)  :: p
      real(xc_f90_kind),    intent(in)  :: rho
      real(xc_f90_kind),    intent(in)  :: sigma
      real(xc_f90_kind),    intent(in)  :: lrho
      real(xc_f90_kind),    intent(in)  :: tau
      real(xc_f90_kind),    intent(out) :: v2rho2
      real(xc_f90_kind),    intent(out) :: v2rhosigma
      real(xc_f90_kind),    intent(out) :: v2sigma2
      real(xc_f90_kind),    intent(out) :: v2rhotau
      real(xc_f90_kind),    intent(out) :: v2tausigma
      real(xc_f90_kind),    intent(out) :: v2tau2
    end subroutine XC_F90(mgga_fxc)
  end interface

  interface
    subroutine XC_F90(mgga_x_tb09_set_par)(p, cc)
      use XC_F90(types_m)
      type(XC_F90(func_t)), intent(inout)  :: p
      real(xc_f90_kind),    intent(in)     :: cc
    end subroutine XC_F90(mgga_x_tb09_set_par)
  end interface


  ! the LCAs
  !----------------------------------------------------------------
  interface
    subroutine XC_F90(lca_init)(p, info, functional, nspin)
      use XC_F90(types_m)
      type(XC_F90(func_t)), intent(out) :: p
      type(XC_F90(info_t)), intent(out) :: info
      integer,         intent(in)  :: functional
      integer,         intent(in)  :: nspin
    end subroutine XC_F90(lca_init)

    subroutine XC_F90(lca_end)(p)
      use XC_F90(types_m)
      type(XC_F90(func_t)), intent(inout) :: p
    end subroutine XC_F90(lca_end)
  end interface


  !----------------------------------------------------------------
  interface
    subroutine XC_F90(lca)(p, rho, v, e, dedd, dedv)
      use XC_F90(types_m)
      type(XC_F90(func_t)),   intent(in)  :: p
      real(xc_f90_kind), intent(in)  :: rho   ! rho(nspin) the density
      real(xc_f90_kind), intent(in)  :: v     ! v(3,nspin) the vorticity
      real(xc_f90_kind), intent(out) :: e     ! the energy per unit particle
      real(xc_f90_kind), intent(out) :: dedd  ! dedd(nspin) the derivative of the energy
                                              ! in terms of the density
      real(xc_f90_kind), intent(out) :: dedv  ! and in terms of the vorticity
    end subroutine XC_F90(lca)
  end interface

end module XC_F90(lib_m)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
