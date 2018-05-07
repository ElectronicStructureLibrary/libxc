!! Copyright (C) 2003-2015 Miguel Marques
!! All rights reserved.
!!
!! This file is dual-licensed under a GPL and a BSD license
!!
!! MPL License:
!!
!! This Source Code Form is subject to the terms of the Mozilla Public
!! License, v. 2.0. If a copy of the MPL was not distributed with this
!! file, You can obtain one at mozilla.org/MPL/2.0/.
!!
!! BSD License:
!!
!! Redistribution and use in source and binary forms, with or without
!! modification, are permitted provided that the following conditions
!! are met:
!!
!! 1. Redistributions of source code must retain the above copyright
!! notice, this list of conditions and the following disclaimer.
!!
!! 2. Redistributions in binary form must reproduce the above
!! copyright notice, this list of conditions and the following
!! disclaimer in the documentation and/or other materials provided
!! with the distribution.
!!
!! 3. Neither the name of the copyright holder nor the names of its
!! contributors may be used to endorse or promote products derived
!! from this software without specific prior written permission.
!!
!! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
!! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
!! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
!! COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
!! INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
!! (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
!! SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
!! HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
!! STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!! ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
!! OF THE POSSIBILITY OF SUCH DAMAGE.


#ifdef SINGLE_PRECISION
#  define xc_f90_x xc_s_f90_ ## x
#else
#  define xc_f90_x xc_f90_ ## x
#endif


!-------------------------------------------------------------------
module xc_f90_types_m
  implicit none
  integer, public, parameter :: xc_f90_kind = selected_real_kind(14)

  type xc_f90_pointer_t
    private
    integer, pointer :: buffer
  end type xc_f90_pointer_t

end module xc_f90_types_m


!-------------------------------------------------------------------
module xc_f90_lib_m

  use xc_f90_types_m
  use libxc_funcs_m

  implicit none

  public

  integer, parameter ::             &
    XC_UNPOLARIZED          =   1,  &  ! Spin unpolarized
    XC_POLARIZED            =   2      ! Spin polarized

  integer, parameter ::             &
    XC_NON_RELATIVISTIC     =   0,  &  ! Functional includes or not relativistic
    XC_RELATIVISTIC         =   1      ! corrections. Only available in some functionals.

  ! Kinds
  integer, parameter ::             &
    XC_EXCHANGE             =   0,  &
    XC_CORRELATION          =   1,  &
    XC_EXCHANGE_CORRELATION =   2,  &
    XC_KINETIC              =   3

  ! Families of xc functionals
  integer, parameter ::     &
    XC_FAMILY_UNKNOWN       =  -1,  &
    XC_FAMILY_NONE          =   0,  &
    XC_FAMILY_LDA           =   1,  &
    XC_FAMILY_GGA           =   2,  &
    XC_FAMILY_MGGA          =   4,  &
    XC_FAMILY_LCA           =   8,  &
    XC_FAMILY_OEP           =  16,  &
    XC_FAMILY_HYB_GGA       =  32,  &
    XC_FAMILY_HYB_MGGA       = 64

  integer, parameter ::                &
    XC_FLAGS_HAVE_EXC        =     1,  &
    XC_FLAGS_HAVE_VXC        =     2,  &
    XC_FLAGS_HAVE_FXC        =     4,  &
    XC_FLAGS_HAVE_KXC        =     8,  &
    XC_FLAGS_HAVE_LXC        =    16,  &
    XC_FLAGS_1D              =    32,  &
    XC_FLAGS_2D              =    64,  &
    XC_FLAGS_3D              =   128,  &
    XC_FLAGS_HYB_CAM         =   256,  &
    XC_FLAGS_HYB_CAMY        =   512,  &
    XC_FLAGS_HYB_VV10        =  1024,  &
    XC_FLAGS_HYB_LC          =  2048,  &
    XC_FLAGS_HYB_LCY         =  4096,  &
    XC_FLAGS_STABLE          =  8192,  &
    XC_FLAGS_DEVELOPMENT     = 16384,  &
    XC_FLAGS_NEEDS_LAPLACIAN = 32768

  integer, parameter, public ::        &
    XC_TAU_EXPLICIT         =     0,   &
    XC_TAU_EXPANSION        =     1

  integer, parameter, public ::        &
    XC_MAX_REFERENCES       =     5

  ! These are old names kept for compatibility
  integer, parameter ::                &
    XC_GGA_X_BGCP           =  38,     &
    XC_GGA_C_BGCP           =  39,     &
    XC_GGA_C_VPBE           =  83,     &
    XC_GGA_XC_LB            = 160,     &
    XC_MGGA_C_CC06          = 229,     &
    XC_GGA_K_ABSR1          = 506,     &
    XC_GGA_K_ABSR2          = 507


  !----------------------------------------------------------------
  interface
    subroutine xc_f90_version(major, minor, micro)
      integer, intent(out) :: major, minor, micro
    end subroutine xc_f90_version

    subroutine xc_f90_version_string(version)
      character(len=10), intent(out) :: version
    end subroutine xc_f90_version_string
  end interface

  !----------------------------------------------------------------
  interface
    integer function xc_f90_info_number(info)
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in) :: info
    end function xc_f90_info_number

    integer function xc_f90_info_kind(info)
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in) :: info
    end function xc_f90_info_kind

    subroutine xc_f90_info_name(info, s)
      use xc_f90_types_m
      type(xc_f90_pointer_t),  intent(in)  :: info
      character(len=*), intent(out) :: s
    end subroutine xc_f90_info_name

    integer function xc_f90_info_family(info)
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in)  :: info
    end function xc_f90_info_family

    integer function xc_f90_info_flags(info)
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in)  :: info
    end function xc_f90_info_flags

    subroutine xc_f90_info_refs(info, number, s)
      use xc_f90_types_m
      type(xc_f90_pointer_t),  intent(in)    :: info
      integer,                  intent(inout) :: number ! number of the reference. Must be 0 in the first call
      character(len=*),         intent(out)   :: s      ! the string that is output
    end subroutine xc_f90_info_refs

    subroutine xc_f90_functional_get_name(func_number, func_string)
      integer, intent(in) :: func_number
      character(len=256), intent(out) :: func_string
    end subroutine xc_f90_functional_get_name

    integer function xc_f90_functional_get_number(func_string)
      character(len=*), intent(in) :: func_string
    end function xc_f90_functional_get_number

    integer function xc_f90_family_from_id(id)
      use xc_f90_types_m
      integer, intent(in) :: id
    end function xc_f90_family_from_id

    integer function xc_f90_number_of_functionals()
      use xc_f90_types_m
    end function xc_f90_number_of_functionals

    integer function xc_f90_maximum_name_length()
      use xc_f90_types_m
    end function xc_f90_maximum_name_length

    subroutine xc_f90_available_functional_numbers(list)
      use xc_f90_types_m
      integer, intent(out) :: list
    end subroutine xc_f90_available_functional_numbers

  end interface


  !----------------------------------------------------------------
  interface
    subroutine xc_f90_func_init(p, info, functional, nspin)
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(out) :: p
      type(xc_f90_pointer_t), intent(out) :: info
      integer,                 intent(in)  :: functional
      integer,                 intent(in)  :: nspin
    end subroutine xc_f90_func_init

    subroutine xc_f90_func_end(p)
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(inout) :: p
    end subroutine xc_f90_func_end

    subroutine xc_f90_func_set_dens_threshold(p, dens_threshold)
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(inout) :: p
      real(xc_f90_kind),       intent(in)    :: dens_threshold
    end subroutine xc_f90_func_set_dens_threshold

    subroutine xc_f90_func_set_ext_params(p, ext_params)
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(inout) :: p
      real(xc_f90_kind),       intent(in)    :: ext_params
    end subroutine xc_f90_func_set_ext_params
  end interface


  ! LDAs
  !----------------------------------------------------------------
  interface
    subroutine xc_f90_lda(p, np, rho, zk, vrho, fxc, kxc)
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in)  :: p
      integer,                 intent(in)  :: np
      real(xc_f90_kind),       intent(in)  :: rho   ! rho(nspin) the density
      real(xc_f90_kind),       intent(out) :: zk    ! the energy per unit particle
      real(xc_f90_kind),       intent(out) :: vrho  ! v(nspin) the potential
      real(xc_f90_kind),       intent(out) :: fxc   ! v(nspin,nspin) the xc kernel
      real(xc_f90_kind),       intent(out) :: kxc   ! v(nspin,nspin,nspin) the derivative of xc kernel
    end subroutine xc_f90_lda

    subroutine xc_f90_lda_exc(p, np, rho, zk)
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in)  :: p
      integer,                 intent(in)  :: np
      real(xc_f90_kind),       intent(in)  :: rho   ! rho(nspin) the density
      real(xc_f90_kind),       intent(out) :: zk    ! the energy per unit particle
    end subroutine xc_f90_lda_exc

    subroutine xc_f90_lda_exc_vxc(p, np, rho, e, v)
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in)  :: p
      integer,                 intent(in)  :: np
      real(xc_f90_kind),       intent(in)  :: rho   ! rho(nspin) the density
      real(xc_f90_kind),       intent(out) :: e     ! the energy per unit particle
      real(xc_f90_kind),       intent(out) :: v     ! v(nspin) the potential
    end subroutine xc_f90_lda_exc_vxc

    subroutine xc_f90_lda_vxc(p, np, rho, v)
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in)  :: p
      integer,                 intent(in)  :: np
      real(xc_f90_kind),       intent(in)  :: rho   ! rho(nspin) the density
      real(xc_f90_kind),       intent(out) :: v     ! v(nspin) the potential
    end subroutine xc_f90_lda_vxc

    subroutine xc_f90_lda_vxc_fxc(p, np, rho, v, fxc)
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in)  :: p
      integer,                 intent(in)  :: np
      real(xc_f90_kind),       intent(in)  :: rho   ! rho(nspin) the density
      real(xc_f90_kind),       intent(out) :: v     ! v(nspin) the potential
      real(xc_f90_kind),       intent(out) :: fxc   ! v(nspin,nspin) the xc kernel
    end subroutine xc_f90_lda_vxc_fxc

    subroutine xc_f90_lda_fxc(p, np, rho, fxc)
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in)  :: p
      integer,                 intent(in)  :: np
      real(xc_f90_kind),       intent(in)  :: rho   ! rho(nspin) the density
      real(xc_f90_kind),       intent(out) :: fxc   ! v(nspin,nspin) the xc kernel
    end subroutine xc_f90_lda_fxc

    subroutine xc_f90_lda_kxc(p, np, rho, kxc)
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in)  :: p
      integer,                 intent(in)  :: np
      real(xc_f90_kind),       intent(in)  :: rho   ! rho(nspin) the density
      real(xc_f90_kind),       intent(out) :: kxc
    end subroutine xc_f90_lda_kxc
  end interface

  ! GGAs
  !----------------------------------------------------------------
  interface
    subroutine xc_f90_gga(p, np, rho, sigma, zk, vrho, vsigma, &
        v2rho2, v2rhosigma, v2sigma2, v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3)
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in)  :: p
      integer,                 intent(in)  :: np
      real(xc_f90_kind),       intent(in)  :: rho
      real(xc_f90_kind),       intent(in)  :: sigma
      real(xc_f90_kind),       intent(out) :: zk
      real(xc_f90_kind),       intent(out) :: vrho
      real(xc_f90_kind),       intent(out) :: vsigma
      real(xc_f90_kind),       intent(out) :: v2rho2
      real(xc_f90_kind),       intent(out) :: v2rhosigma
      real(xc_f90_kind),       intent(out) :: v2sigma2
      real(xc_f90_kind),       intent(out) :: v3rho3
      real(xc_f90_kind),       intent(out) :: v3rho2sigma
      real(xc_f90_kind),       intent(out) :: v3rhosigma2
      real(xc_f90_kind),       intent(out) :: v3sigma3
    end subroutine xc_f90_gga

    subroutine xc_f90_gga_exc(p, np, rho, sigma, zk)
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in)  :: p
      integer,                 intent(in)  :: np
      real(xc_f90_kind),       intent(in)  :: rho
      real(xc_f90_kind),       intent(in)  :: sigma
      real(xc_f90_kind),       intent(out) :: zk
    end subroutine xc_f90_gga_exc

    subroutine xc_f90_gga_exc_vxc(p, np, rho, sigma, zk, vrho, vsigma)
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in)  :: p
      integer,                 intent(in)  :: np
      real(xc_f90_kind),       intent(in)  :: rho
      real(xc_f90_kind),       intent(in)  :: sigma
      real(xc_f90_kind),       intent(out) :: zk
      real(xc_f90_kind),       intent(out) :: vrho
      real(xc_f90_kind),       intent(out) :: vsigma
    end subroutine xc_f90_gga_exc_vxc

    subroutine xc_f90_gga_vxc(p, np, rho, sigma, vrho, vsigma)
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in)  :: p
      integer,                 intent(in)  :: np
      real(xc_f90_kind),       intent(in)  :: rho
      real(xc_f90_kind),       intent(in)  :: sigma
      real(xc_f90_kind),       intent(out) :: vrho
      real(xc_f90_kind),       intent(out) :: vsigma
    end subroutine xc_f90_gga_vxc

    subroutine xc_f90_gga_vxc_fxc(p, np, rho, sigma, vrho, vsigma, v2rho2, v2rhosigma, v2sigma2)
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in)  :: p
      integer,                 intent(in)  :: np
      real(xc_f90_kind),       intent(in)  :: rho
      real(xc_f90_kind),       intent(in)  :: sigma
      real(xc_f90_kind),       intent(out) :: vrho
      real(xc_f90_kind),       intent(out) :: vsigma
      real(xc_f90_kind),       intent(out) :: v2rho2
      real(xc_f90_kind),       intent(out) :: v2rhosigma
      real(xc_f90_kind),       intent(out) :: v2sigma2
    end subroutine xc_f90_gga_vxc_fxc

    subroutine xc_f90_gga_fxc(p, np, rho, sigma, v2rho2, v2rhosigma, v2sigma2)
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in)  :: p
      integer,                 intent(in)  :: np
      real(xc_f90_kind),       intent(in)  :: rho
      real(xc_f90_kind),       intent(in)  :: sigma
      real(xc_f90_kind),       intent(out) :: v2rho2
      real(xc_f90_kind),       intent(out) :: v2rhosigma
      real(xc_f90_kind),       intent(out) :: v2sigma2
    end subroutine xc_f90_gga_fxc

    subroutine xc_f90_gga_kxc(p, np, rho, sigma, v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3)
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in)  :: p
      integer,                 intent(in)  :: np
      real(xc_f90_kind),       intent(in)  :: rho
      real(xc_f90_kind),       intent(in)  :: sigma
      real(xc_f90_kind),       intent(out) :: v3rho3
      real(xc_f90_kind),       intent(out) :: v3rho2sigma
      real(xc_f90_kind),       intent(out) :: v3rhosigma2
      real(xc_f90_kind),       intent(out) :: v3sigma3
    end subroutine xc_f90_gga_kxc
  end interface


  !----------------------------------------------------------------
  interface
    subroutine xc_f90_gga_lb_modified(p, np, rho, grho, r, dedd)
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in)  :: p
      integer,              intent(in)  :: np
      real(xc_f90_kind),    intent(in)  :: rho   ! rho(nspin) the density
      real(xc_f90_kind),    intent(in)  :: grho  ! grho(3,nspin) the gradient of the density
      real(xc_f90_kind),    intent(in)  :: r     ! distance from center of finite system
      real(xc_f90_kind),    intent(out) :: dedd
    end subroutine xc_f90_gga_lb_modified
  end interface


  !----------------------------------------------------------------
  interface
    subroutine xc_f90_gga_ak13_get_asymptotic(homo, asymp)
      use xc_f90_types_m
      real(xc_f90_kind),       intent(in)  :: homo
      real(xc_f90_kind),       intent(out) :: asymp
    end subroutine xc_f90_gga_ak13_get_asymptotic
  end interface

  !----------------------------------------------------------------
  interface
    subroutine xc_f90_hyb_exx_coef(p, coef)
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in)  :: p
      real(xc_f90_kind),       intent(out) :: coef
    end subroutine xc_f90_hyb_exx_coef

    subroutine xc_f90_hyb_cam_coef(p, omega, alpha, beta)
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in)  :: p
      real(xc_f90_kind),       intent(out) :: omega, alpha, beta
    end subroutine xc_f90_hyb_cam_coef

    subroutine xc_f90_nlc_coef(p, nlc_b, nlc_c)
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in)  :: p
      real(xc_f90_kind),       intent(out) :: nlc_b, nlc_c
    end subroutine xc_f90_nlc_coef
  end interface


  ! the meta-GGAs
  !----------------------------------------------------------------
  interface
    subroutine xc_f90_mgga(p, np, rho, sigma, lapl, tau, zk, vrho, vsigma, vlapl, vtau, &
      v2rho2, v2sigma2, v2lapl2, v2tau2, v2rhosigma, v2rholapl, v2rhotau, &
      v2sigmalapl, v2sigmatau, v2lapltau)

      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in)  :: p
      integer,              intent(in)  :: np
      real(xc_f90_kind),    intent(in)  :: rho
      real(xc_f90_kind),    intent(in)  :: sigma
      real(xc_f90_kind),    intent(in)  :: lapl
      real(xc_f90_kind),    intent(in)  :: tau
      real(xc_f90_kind),    intent(out) :: zk
      real(xc_f90_kind),    intent(out) :: vrho
      real(xc_f90_kind),    intent(out) :: vsigma
      real(xc_f90_kind),    intent(out) :: vlapl
      real(xc_f90_kind),    intent(out) :: vtau
      real(xc_f90_kind),    intent(out) :: v2rho2
      real(xc_f90_kind),    intent(out) :: v2sigma2
      real(xc_f90_kind),    intent(out) :: v2lapl2
      real(xc_f90_kind),    intent(out) :: v2tau2
      real(xc_f90_kind),    intent(out) :: v2rhosigma
      real(xc_f90_kind),    intent(out) :: v2rholapl
      real(xc_f90_kind),    intent(out) :: v2rhotau
      real(xc_f90_kind),    intent(out) :: v2sigmalapl
      real(xc_f90_kind),    intent(out) :: v2sigmatau
      real(xc_f90_kind),    intent(out) :: v2lapltau
    end subroutine xc_f90_mgga

    subroutine xc_f90_mgga_exc(p, np, rho, sigma, lapl, tau, zk)
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in)  :: p
      integer,              intent(in)  :: np
      real(xc_f90_kind),    intent(in)  :: rho
      real(xc_f90_kind),    intent(in)  :: sigma
      real(xc_f90_kind),    intent(in)  :: lapl
      real(xc_f90_kind),    intent(in)  :: tau
      real(xc_f90_kind),    intent(out) :: zk
    end subroutine xc_f90_mgga_exc

    subroutine xc_f90_mgga_exc_vxc(p, np, rho, sigma, lapl, tau, zk, vrho, vsigma, vlapl, vtau)
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in)  :: p
      integer,              intent(in)  :: np
      real(xc_f90_kind),    intent(in)  :: rho
      real(xc_f90_kind),    intent(in)  :: sigma
      real(xc_f90_kind),    intent(in)  :: lapl
      real(xc_f90_kind),    intent(in)  :: tau
      real(xc_f90_kind),    intent(out) :: zk
      real(xc_f90_kind),    intent(out) :: vrho
      real(xc_f90_kind),    intent(out) :: vsigma
      real(xc_f90_kind),    intent(out) :: vlapl
      real(xc_f90_kind),    intent(out) :: vtau
    end subroutine xc_f90_mgga_exc_vxc

    subroutine xc_f90_mgga_vxc(p, np, rho, sigma, lapl, tau, vrho, vsigma, vlapl, vtau)
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in)  :: p
      integer,              intent(in)  :: np
      real(xc_f90_kind),    intent(in)  :: rho
      real(xc_f90_kind),    intent(in)  :: sigma
      real(xc_f90_kind),    intent(in)  :: lapl
      real(xc_f90_kind),    intent(in)  :: tau
      real(xc_f90_kind),    intent(out) :: vrho
      real(xc_f90_kind),    intent(out) :: vsigma
      real(xc_f90_kind),    intent(out) :: vlapl
      real(xc_f90_kind),    intent(out) :: vtau
    end subroutine xc_f90_mgga_vxc

    subroutine xc_f90_mgga_vxc_fxc(p, np, rho, sigma, lapl, tau, &
      vrho, vsigma, vlapl, vtau, &
      v2rho2, v2sigma2, v2lapl2, v2tau2, v2rhosigma, v2rholapl, v2rhotau, &
      v2sigmalapl, v2sigmatau, v2lapltau)

      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in)  :: p
      integer,              intent(in)  :: np
      real(xc_f90_kind),    intent(in)  :: rho
      real(xc_f90_kind),    intent(in)  :: sigma
      real(xc_f90_kind),    intent(in)  :: lapl
      real(xc_f90_kind),    intent(in)  :: tau
      real(xc_f90_kind),    intent(out) :: vrho
      real(xc_f90_kind),    intent(out) :: vsigma
      real(xc_f90_kind),    intent(out) :: vlapl
      real(xc_f90_kind),    intent(out) :: vtau
      real(xc_f90_kind),    intent(out) :: v2rho2
      real(xc_f90_kind),    intent(out) :: v2sigma2
      real(xc_f90_kind),    intent(out) :: v2lapl2
      real(xc_f90_kind),    intent(out) :: v2tau2
      real(xc_f90_kind),    intent(out) :: v2rhosigma
      real(xc_f90_kind),    intent(out) :: v2rholapl
      real(xc_f90_kind),    intent(out) :: v2rhotau
      real(xc_f90_kind),    intent(out) :: v2sigmalapl
      real(xc_f90_kind),    intent(out) :: v2sigmatau
      real(xc_f90_kind),    intent(out) :: v2lapltau
    end subroutine xc_f90_mgga_vxc_fxc

    subroutine xc_f90_mgga_fxc(p, np, rho, sigma, lapl, tau, &
      v2rho2, v2sigma2, v2lapl2, v2tau2, v2rhosigma, v2rholapl, v2rhotau, &
      v2sigmalapl, v2sigmatau, v2lapltau)

      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in)  :: p
      integer,              intent(in)  :: np
      real(xc_f90_kind),    intent(in)  :: rho
      real(xc_f90_kind),    intent(in)  :: sigma
      real(xc_f90_kind),    intent(in)  :: lapl
      real(xc_f90_kind),    intent(in)  :: tau
      real(xc_f90_kind),    intent(out) :: v2rho2
      real(xc_f90_kind),    intent(out) :: v2sigma2
      real(xc_f90_kind),    intent(out) :: v2lapl2
      real(xc_f90_kind),    intent(out) :: v2tau2
      real(xc_f90_kind),    intent(out) :: v2rhosigma
      real(xc_f90_kind),    intent(out) :: v2rholapl
      real(xc_f90_kind),    intent(out) :: v2rhotau
      real(xc_f90_kind),    intent(out) :: v2sigmalapl
      real(xc_f90_kind),    intent(out) :: v2sigmatau
      real(xc_f90_kind),    intent(out) :: v2lapltau
    end subroutine xc_f90_mgga_fxc
  end interface

end module xc_f90_lib_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
