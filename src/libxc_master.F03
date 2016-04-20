!! Copyright (C) 2016 Micael Oliveira
!! All rights reserved.
!!
!! This file is dual-licensed under a GPL and a BSD license
!!
!! GPL License:
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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
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

!! $Id$


#ifdef SINGLE_PRECISION
#  define XC(x) xc_s_ ## x
#  define FXC(x) xc_s_f03_ ## x
#  define RTYPE c_float
#else
#  define XC(x) xc_ ## x
#  define FXC(x) xc_f03_ ## x
#  define RTYPE c_double
#endif

module FXC(lib_m)
  use, intrinsic :: iso_c_binding
  implicit none

  private
  public :: &
    ! version
    FXC(version), &
    FXC(version_string), &
    ! info
    FXC(func_info_t), &
    FXC(func_info_get_number), &
    FXC(func_info_get_kind), &
    FXC(func_info_get_name), &
    FXC(func_info_get_family), &
    FXC(func_info_get_flags), &
    FXC(func_info_get_refs), &
    ! func
    FXC(func_t), &
    FXC(func_init), &
    FXC(func_end), &
    FXC(func_get_info), &
    FXC(functional_get_name), &
    FXC(functional_get_number), &
    FXC(family_from_id), &
    ! lda
    FXC(lda), &
    FXC(lda_exc), &
    FXC(lda_exc_vxc), &
    FXC(lda_vxc), &
    FXC(lda_fxc), &
    FXC(lda_kxc), &
    FXC(lda_x_set_params), &
    FXC(lda_x_1d_set_params), &
    FXC(lda_c_1d_csc_set_params), &
    FXC(lda_c_xalpha_set_params), &
    FXC(lda_c_2d_prm_set_params), &
    FXC(lda_c_vwn_set_params), &
    FXC(lda_xc_ksdt_set_params), &
    ! gga
    FXC(gga), &
    FXC(gga_exc), &
    FXC(gga_exc_vxc), &
    FXC(gga_vxc), &
    FXC(gga_fxc), &
    FXC(gga_kxc), &
    FXC(gga_lb_modified), &
    FXC(gga_x_b86_set_params), &
    FXC(gga_x_b88_set_params), &
    FXC(gga_x_pbe_set_params), &
    FXC(gga_c_pbe_set_params), &
    FXC(gga_x_pw91_set_params), &
    FXC(gga_x_pw91_set_params2), &
    FXC(gga_x_rpbe_set_params), &
    FXC(gga_x_optx_set_params), &
    FXC(gga_c_lyp_set_params), &
    FXC(gga_lb_set_params), &
    FXC(gga_k_tflw_set_params), &
    FXC(gga_x_2d_b88_set_params), &
    FXC(gga_x_wpbeh_set_params), &
    FXC(gga_x_hjs_set_params), &
    FXC(gga_x_ityh_set_params), &
    FXC(gga_x_sfat_set_params), &
    FXC(gga_x_ssb_sw_set_params), &
    FXC(gga_x_kt_set_params), &
    FXC(gga_x_lambda_set_params), &
    FXC(gga_ak13_get_asymptotic), &
    FXC(hyb_exx_coef), &
    FXC(hyb_cam_coef), &
    FXC(nlc_coef), &
    FXC(hyb_gga_xc_hse_set_params), &
    FXC(hyb_gga_xc_pbeh_set_params), &
    ! mgga
    FXC(mgga), &
    FXC(mgga_exc), &
    FXC(mgga_exc_vxc), &
    FXC(mgga_vxc), &
    FXC(mgga_fxc), &
    FXC(mgga_x_tb09_set_params), &
    FXC(mgga_x_tpss_set_params), &
    FXC(mgga_c_bc95_set_params), &
    FXC(mgga_c_pkzb_set_params)


  integer(c_int), parameter, public :: &
    XC_UNPOLARIZED          =   1,     &  ! Spin unpolarized
    XC_POLARIZED            =   2         ! Spin polarized

    integer(c_int), parameter, public :: &
    XC_NON_RELATIVISTIC     =   0,     &  ! Functional includes or not relativistic
    XC_RELATIVISTIC         =   1         ! corrections. Only available in some functionals.

  ! Kinds
  integer(c_int), parameter, public :: &
    XC_EXCHANGE             =   0,     &
    XC_CORRELATION          =   1,     &
    XC_EXCHANGE_CORRELATION =   2,     &
    XC_KINETIC              =   3

  ! Families of xc functionals
  integer(c_int), parameter, public :: &
    XC_FAMILY_UNKNOWN       =  -1,     &
    XC_FAMILY_NONE          =   0,     &
    XC_FAMILY_LDA           =   1,     &
    XC_FAMILY_GGA           =   2,     &
    XC_FAMILY_MGGA          =   4,     &
    XC_FAMILY_LCA           =   8,     &
    XC_FAMILY_OEP           =  16,     &
    XC_FAMILY_HYB_GGA       =  32,     &
    XC_FAMILY_HYB_MGGA      =  64

  integer(c_int), parameter, public :: &
    XC_FLAGS_HAVE_EXC       =     1,   &
    XC_FLAGS_HAVE_VXC       =     2,   &
    XC_FLAGS_HAVE_FXC       =     4,   &
    XC_FLAGS_HAVE_KXC       =     8,   &
    XC_FLAGS_HAVE_LXC       =    16,   &
    XC_FLAGS_1D             =    32,   &
    XC_FLAGS_2D             =    64,   &
    XC_FLAGS_3D             =   128,   &
    XC_FLAGS_HYB_CAM        =   256,   &
    XC_FLAGS_HYB_CAMY       =   512,   &
    XC_FLAGS_VV10           =  1024,   &
    XC_FLAGS_HYB_LC         =  2048,   &
    XC_FLAGS_HYB_LCY        =  4096,   &
    XC_FLAGS_STABLE         =  8192,   &
    XC_FLAGS_DEVELOPMENT    = 16384

  integer(c_int), parameter, public :: &
    XC_TAU_EXPLICIT         =     0,   &
    XC_TAU_EXPANSION        =     1

  ! List of functionals
#include "libxc_inc.f03"
  
  ! These are old names keep for compatibility, and that should disappear soon
  integer(c_int), parameter, public :: &
    XC_GGA_C_VPBE           =  83,     &
    XC_GGA_XC_LB            = 160,     &
    XC_GGA_K_ABSR1          = 506,     &
    XC_GGA_K_ABSR2          = 507

  !----------------------------------------------------------------
  interface
    subroutine XC(version)(major, minor, micro) bind(c)
      import
      integer(c_int), intent(out) :: major, minor, micro
    end subroutine XC(version)

    type(c_ptr) function XC(version_string)() bind(c)
      import
    end function XC(version_string)
  end interface

  
  !----------------------------------------------------------------
  type :: FXC(func_info_t)
    private
    type(c_ptr) :: ptr = C_NULL_PTR
  end type FXC(func_info_t)

  interface
    integer(c_int) function XC(func_info_get_number)(info) bind(c)
      import
      type(c_ptr), value :: info
    end function XC(func_info_get_number)
    
    integer(c_int) function XC(func_info_get_kind)(info) bind(c)
      import
      type(c_ptr), value :: info
    end function XC(func_info_get_kind)

    type(c_ptr) function XC(func_info_get_name)(info) bind(c)
      import
      type(c_ptr), value :: info
    end function XC(func_info_get_name)

    integer(c_int) function XC(func_info_get_family)(info) bind(c)
      import
      type(c_ptr), value :: info
    end function XC(func_info_get_family)

    integer(c_int) function XC(func_info_get_flags)(info) bind(c)
      import
      type(c_ptr), value :: info
    end function XC(func_info_get_flags)

    type(c_ptr) function XC(func_info_get_ref)(info, number) bind(c)
      import
      type(c_ptr),    value :: info
      integer(c_int), value :: number
    end function XC(func_info_get_ref)

  end interface


  !----------------------------------------------------------------
  type :: FXC(func_t)
    private
    type(c_ptr) :: ptr = C_NULL_PTR
  end type FXC(func_t)

  interface
    type(c_ptr) function XC(func_alloc)() bind(c)
      import
    end function XC(func_alloc)
    
    integer(c_int) function XC(func_init)(p, functional, nspin) bind(c)
      import
      type(c_ptr),    value :: p
      integer(c_int), value :: functional, nspin
    end function XC(func_init)

    subroutine XC(func_end)(p) bind(c)
      import
      type(c_ptr), value :: p
    end subroutine XC(func_end)

    subroutine XC(func_free)(p) bind(c)
      import
      type(c_ptr), value :: p
    end subroutine XC(func_free)

    type(c_ptr) function XC(func_get_info)(p) bind(c)
      import
      type(c_ptr), value :: p
    end function XC(func_get_info)
    
    type(c_ptr) function XC(functional_get_name)(number) bind(c)
      import
      integer(c_int), value :: number
    end function XC(functional_get_name)

    integer(c_int) function XC(functional_get_number)(func_string) bind(c)
      import
      character(kind=c_char) :: func_string(*)
    end function XC(functional_get_number)
 
    integer(c_int) function XC(family_from_id)(id, family, number) bind(c)
      import
      integer(c_int), value :: id
      type(c_ptr),    value :: family, number
    end function XC(family_from_id)
  end interface


  ! LDAs
  !----------------------------------------------------------------
  interface
    subroutine XC(lda)(p, np, rho, zk, vrho, v2rho2, v3rho3) bind(c)
      import
      type(c_ptr),    value       :: p
      integer(c_int), value       :: np
      real(RTYPE),    intent(in)  :: rho(*)
      real(RTYPE),    intent(out) :: zk(*), vrho(*), v2rho2(*), v3rho3(*)
    end subroutine XC(lda)

    subroutine XC(lda_exc)(p, np, rho, zk) bind(c)
      import
      type(c_ptr),    value       :: p
      integer(c_int), value       :: np
      real(RTYPE),    intent(in)  :: rho(*)
      real(RTYPE),    intent(out) :: zk(*)
    end subroutine XC(lda_exc)

    subroutine XC(lda_exc_vxc)(p, np, rho, zk, vrho) bind(c)
      import
      type(c_ptr),    value       :: p
      integer(c_int), value       :: np
      real(RTYPE),    intent(in)  :: rho(*)
      real(RTYPE),    intent(out) :: zk(*), vrho(*)
    end subroutine XC(lda_exc_vxc)

    subroutine XC(lda_vxc)(p, np, rho, vrho) bind(c)
      import
      type(c_ptr),    value       :: p
      integer(c_int), value       :: np
      real(RTYPE),    intent(in)  :: rho(*)
      real(RTYPE),    intent(out) :: vrho(*)
    end subroutine XC(lda_vxc)

    subroutine XC(lda_fxc)(p, np, rho, v2rho2) bind(c)
      import
      type(c_ptr),    value       :: p
      integer(c_int), value       :: np
      real(RTYPE),    intent(in)  :: rho(*)
      real(RTYPE),    intent(out) :: v2rho2(*)
    end subroutine XC(lda_fxc)

    subroutine XC(lda_kxc)(p, np, rho, v3rho3) bind(c)
      import
      type(c_ptr),    value       :: p
      integer(c_int), value       :: np
      real(RTYPE),    intent(in)  :: rho(*)
      real(RTYPE),    intent(out) :: v3rho3(*)
    end subroutine XC(lda_kxc)
  end interface
  

  interface
    subroutine XC(lda_x_set_params)(p, alpha, relativistic, omega) bind(c)
      import
      type(c_ptr),    value :: p
      real(RTYPE),    value :: alpha, omega
      integer(c_int), value :: relativistic
    end subroutine XC(lda_x_set_params)

    subroutine XC(lda_x_1d_set_params)(p, interaction, bb) bind(c)
      import
      type(c_ptr),    value  :: p
      integer(c_int), value  :: interaction
      real(RTYPE),    value  :: bb
    end subroutine XC(lda_x_1d_set_params)

    subroutine XC(lda_c_1d_csc_set_params)(p, interaction, bb) bind(c)
      import
      type(c_ptr),    value :: p
      integer(c_int), value :: interaction
      real(RTYPE),    value :: bb
    end subroutine XC(lda_c_1d_csc_set_params)
    
    subroutine XC(lda_c_xalpha_set_params)(p, alpha) bind(c)
      import
      type(c_ptr), value :: p
      real(RTYPE), value :: alpha
    end subroutine XC(lda_c_xalpha_set_params)

    subroutine XC(lda_c_2d_prm_set_params)(p, N) bind(c)
      import
      type(c_ptr), value :: p
      real(RTYPE), value :: N
    end subroutine XC(lda_c_2d_prm_set_params)

    subroutine XC(lda_c_vwn_set_params)(p, spin_interpolation) bind(c)
      import
      type(c_ptr),    value :: p
      integer(c_int), value :: spin_interpolation
    end subroutine XC(lda_c_vwn_set_params)

    subroutine XC(lda_xc_ksdt_set_params)(p, t) bind(c)
      import
      type(c_ptr), value :: p
      real(RTYPE), value :: t
    end subroutine XC(lda_xc_ksdt_set_params)
  end interface

  ! GGAs
  !----------------------------------------------------------------
  interface
    subroutine XC(gga)(p, np, rho, sigma, zk, vrho, vsigma, &
      v2rho2, v2rhosigma, v2sigma2, v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3) bind(c)
      import
      type(c_ptr),    value       :: p
      integer(c_int), value       :: np
      real(RTYPE),    intent(in)  :: rho(*), sigma(*)
      real(RTYPE),    intent(out) :: zk(*), vrho(*), vsigma(*)
      real(RTYPE),    intent(out) :: v2rho2(*), v2rhosigma(*), v2sigma2(*)
      real(RTYPE),    intent(out) :: v3rho3(*), v3rho2sigma(*), v3rhosigma2(*), v3sigma3(*)
    end subroutine XC(gga)

    subroutine XC(gga_exc)(p, np, rho, sigma, zk) bind(c)
      import
      type(c_ptr),    value       :: p
      integer(c_int), value       :: np
      real(RTYPE),    intent(in)  :: rho(*), sigma(*)
      real(RTYPE),    intent(out) :: zk(*)
    end subroutine XC(gga_exc)

    subroutine XC(gga_exc_vxc)(p, np, rho, sigma, zk, vrho, vsigma) bind(c)
      import
      type(c_ptr),    value       :: p
      integer(c_int), value       :: np
      real(RTYPE),    intent(in)  :: rho(*), sigma(*)
      real(RTYPE),    intent(out) :: zk(*), vrho(*), vsigma(*)      
    end subroutine XC(gga_exc_vxc)

    subroutine XC(gga_vxc)(p, np, rho, sigma, vrho, vsigma) bind(c)
      import
      type(c_ptr),    value       :: p
      integer(c_int), value       :: np
      real(RTYPE),    intent(in)  :: rho(*), sigma(*)
      real(RTYPE),    intent(out) :: vrho(*), vsigma(*)
    end subroutine XC(gga_vxc)

    subroutine XC(gga_fxc)(p, np, rho, sigma, v2rho2, v2rhosigma, v2sigma2) bind(c)
      import
      type(c_ptr),    value       :: p
      integer(c_int), value       :: np
      real(RTYPE),    intent(in)  :: rho(*), sigma(*)
      real(RTYPE),    intent(out) :: v2rho2(*), v2rhosigma(*), v2sigma2(*)
    end subroutine XC(gga_fxc)

    subroutine XC(gga_kxc)(p, np, rho, sigma, v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3) bind(c)
      import
      type(c_ptr),    value       :: p
      integer(c_int), value       :: np
      real(RTYPE),    intent(in)  :: rho(*), sigma(*)
      real(RTYPE),    intent(out) :: v3rho3(*), v3rho2sigma(*), v3rhosigma2(*), v3sigma3(*)      
    end subroutine XC(gga_kxc)
  end interface


  interface
    subroutine XC(gga_lb_modified)(p, np, rho, grho, r, dedd) bind(c)
      import
      type(c_ptr),    value       :: p
      integer(c_int), value       :: np
      real(RTYPE),    intent(in)  :: rho(*), grho(*)
      real(RTYPE),    value       :: r
      real(RTYPE),    intent(out) :: dedd(*)
    end subroutine XC(gga_lb_modified)

    subroutine XC(gga_x_b86_set_params)(p, beta, gamma, omega) bind(c)
      import
      type(c_ptr), value :: p
      real(RTYPE), value :: beta, gamma, omega
    end subroutine XC(gga_x_b86_set_params)

    subroutine XC(gga_x_b88_set_params)(p, beta, gamma) bind(c)
      import
      type(c_ptr), value :: p
      real(RTYPE), value :: beta, gamma
    end subroutine XC(gga_x_b88_set_params)

    subroutine XC(gga_x_pbe_set_params)(p, kappa, mu) bind(c)
      import
      type(c_ptr), value :: p
      real(RTYPE), value :: kappa, mu
    end subroutine XC(gga_x_pbe_set_params)

    subroutine XC(gga_c_pbe_set_params)(p, beta) bind(c)
      import
      type(c_ptr), value :: p
      real(RTYPE), value :: beta
    end subroutine XC(gga_c_pbe_set_params)

    subroutine XC(gga_x_pw91_set_params)(p, a, b, c, d, f, alpha, expo) bind(c)
      import
      type(c_ptr), value :: p
      real(RTYPE), value :: a, b, c, d, f, alpha, expo
    end subroutine XC(gga_x_pw91_set_params)

    subroutine XC(gga_x_pw91_set_params2)(p, bt, alpha, expo) bind(c)
      import
      type(c_ptr), value :: p
      real(RTYPE), value :: bt, alpha, expo
    end subroutine XC(gga_x_pw91_set_params2)

    subroutine XC(gga_x_rpbe_set_params)(p, kappa, mu) bind(c)
      import
      type(c_ptr), value :: p
      real(RTYPE), value :: kappa, mu
    end subroutine XC(gga_x_rpbe_set_params)

    subroutine XC(gga_x_optx_set_params)(p, a, b, gamma) bind(c)
      import
      type(c_ptr), value :: p
      real(RTYPE), value :: a, b, gamma
    end subroutine XC(gga_x_optx_set_params)

    subroutine XC(gga_c_lyp_set_params)(p, a, b, c, d) bind(c)
      import
      type(c_ptr), value :: p
      real(RTYPE), value :: a, b, c, d
    end subroutine XC(gga_c_lyp_set_params)

    subroutine XC(gga_lb_set_params)(p, modified, threshold, ip, qtot) bind(c)
      import
      type(c_ptr),    value :: p
      integer(c_int), value :: modified
      real(RTYPE),    value :: threshold, ip, qtot
    end subroutine XC(gga_lb_set_params)

    subroutine XC(gga_k_tflw_set_params)(p, gamma, lambda, n) bind(c)
      import
      type(c_ptr),    value :: p
      real(RTYPE),    value :: gamma, lambda
      integer(c_int), value :: n
    end subroutine XC(gga_k_tflw_set_params)
    
    subroutine XC(gga_x_2d_b88_set_params)(p, beta) bind(c)
      import
      type(c_ptr), value :: p
      real(RTYPE), value :: beta
    end subroutine XC(gga_x_2d_b88_set_params)

    subroutine XC(gga_x_wpbeh_set_params)(p, omega) bind(c)
      import
      type(c_ptr), value :: p
      real(RTYPE), value :: omega
    end subroutine XC(gga_x_wpbeh_set_params)

    subroutine XC(gga_x_hjs_set_params)(p, omega) bind(c)
      import
      type(c_ptr), value :: p
      real(RTYPE), value :: omega
    end subroutine XC(gga_x_hjs_set_params)

    subroutine XC(gga_x_ityh_set_params)(p, func_id, omega) bind(c)
      import
      type(c_ptr),    value :: p
      integer(c_int), value :: func_id
      real(RTYPE),    value :: omega
    end subroutine XC(gga_x_ityh_set_params)

    subroutine XC(gga_x_sfat_set_params)(p, func_id, omega) bind(c)
      import
      type(c_ptr),    value :: p
      integer(c_int), value :: func_id
      real(RTYPE),    value :: omega
    end subroutine XC(gga_x_sfat_set_params)
    
    subroutine XC(gga_x_ssb_sw_set_params)(p, a, b, c, d, e) bind(c)
      import
      type(c_ptr), value :: p
      real(RTYPE), value :: a, b, c, d, e
    end subroutine XC(gga_x_ssb_sw_set_params)

    subroutine XC(gga_x_kt_set_params)(p, gamma, delta) bind(c)
      import
      type(c_ptr), value :: p
      real(RTYPE), value :: gamma, delta
    end subroutine XC(gga_x_kt_set_params)

    subroutine XC(gga_x_lambda_set_params)(p, n) bind(c)
      import
      type(c_ptr),    value :: p
      integer(c_int), value :: n
    end subroutine XC(gga_x_lambda_set_params)

    real(RTYPE) function XC(gga_ak13_get_asymptotic)(homo) bind(c)
      import
      real(RTYPE), value :: homo
    end function XC(gga_ak13_get_asymptotic)
  end interface


  interface
    real(RTYPE) function XC(hyb_exx_coef)(p) bind(c)
      import
      type(c_ptr), value :: p
    end function XC(hyb_exx_coef)

    subroutine XC(hyb_cam_coef)(p, omega, alpha, beta) bind(c)
      import
      type(c_ptr), value       :: p
      real(RTYPE), intent(out) :: omega, alpha, beta
    end subroutine XC(hyb_cam_coef)

    subroutine XC(nlc_coef)(p, nlc_b, nlc_c) bind(c)
      import
      type(c_ptr), value       :: p
      real(RTYPE), intent(out) :: nlc_b, nlc_c
    end subroutine XC(nlc_coef)
    
    subroutine XC(hyb_gga_xc_hse_set_params)(p, alpha, omega) bind(c)
      import
      type(c_ptr), value :: p
      real(RTYPE), value :: alpha, omega
    end subroutine XC(hyb_gga_xc_hse_set_params)

    subroutine XC(hyb_gga_xc_pbeh_set_params)(p, alpha) bind(c)
      import
      type(c_ptr), value :: p
      real(RTYPE), value :: alpha
    end subroutine XC(hyb_gga_xc_pbeh_set_params)
  end interface


  ! the meta-GGAs
  !----------------------------------------------------------------
  interface
    subroutine XC(mgga)(p, np, rho, sigma, lapl, tau, zk, vrho, vsigma, vlapl, vtau, &
      v2rho2, v2sigma2, v2lapl2, v2tau2, v2rhosigma, v2rholapl, v2rhotau, &
      v2sigmalapl, v2sigmatau, v2lapltau) bind(c)
      import
      type(c_ptr),    value       :: p
      integer(c_int), value       :: np
      real(RTYPE),    intent(in)  :: rho(*), sigma(*), lapl(*), tau(*)
      real(RTYPE),    intent(out) :: zk(*), vrho(*), vsigma(*), vlapl(*), vtau(*)
      real(RTYPE),    intent(out) :: v2rho2(*), v2sigma2(*), v2lapl2(*), v2tau2(*), v2rhosigma(*), v2rholapl(*), &
                                     v2rhotau(*), v2sigmalapl(*), v2sigmatau(*), v2lapltau(*)
    end subroutine XC(mgga)

    subroutine XC(mgga_exc)(p, np, rho, sigma, lapl, tau, zk) bind(c)
      import
      type(c_ptr),    value       :: p
      integer(c_int), value       :: np
      real(RTYPE),    intent(in)  :: rho(*), sigma(*), lapl(*), tau(*)
      real(RTYPE),    intent(out) :: zk(*)
    end subroutine XC(mgga_exc)

    subroutine XC(mgga_exc_vxc)(p, np, rho, sigma, lapl, tau, zk, vrho, vsigma, vlapl, vtau) bind(c)
      import
      type(c_ptr),    value       :: p
      integer(c_int), value       :: np
      real(RTYPE),    intent(in)  :: rho(*), sigma(*), lapl(*), tau(*)
      real(RTYPE),    intent(out) :: zk(*), vrho(*), vsigma(*), vlapl(*), vtau(*)
    end subroutine XC(mgga_exc_vxc)

    subroutine XC(mgga_vxc)(p, np, rho, sigma, lapl, tau, vrho, vsigma, vlapl, vtau) bind(c)
      import
      type(c_ptr),    value       :: p
      integer(c_int), value       :: np
      real(RTYPE),    intent(in)  :: rho(*), sigma(*), lapl(*), tau(*)
      real(RTYPE),    intent(out) :: vrho(*), vsigma(*), vlapl(*), vtau(*)
    end subroutine XC(mgga_vxc)

    subroutine XC(mgga_fxc)(p, np, rho, sigma, lapl, tau, &
      v2rho2, v2sigma2, v2lapl2, v2tau2, v2rhosigma, v2rholapl, v2rhotau, &
      v2sigmalapl, v2sigmatau, v2lapltau) bind(c)
      import
      type(c_ptr),    value       :: p
      integer(c_int), value       :: np
      real(RTYPE),    intent(in)  :: rho(*), sigma(*), lapl(*), tau(*)
      real(RTYPE),    intent(out) :: v2rho2(*), v2sigma2(*), v2lapl2(*), v2tau2(*), v2rhosigma(*), v2rholapl(*), &
                                     v2rhotau(*), v2sigmalapl(*), v2sigmatau(*), v2lapltau(*)
    end subroutine XC(mgga_fxc)
  end interface

  interface
    subroutine XC(mgga_x_tb09_set_params)(p, c) bind(c)
      import
      type(c_ptr), value :: p
      real(RTYPE), value :: c
    end subroutine XC(mgga_x_tb09_set_params)

    subroutine XC(mgga_x_tpss_set_params)(p, b, c, e, kappa, mu) bind(c)
      import
      type(c_ptr), value :: p
      real(RTYPE), value :: b, c, e, kappa, mu
    end subroutine XC(mgga_x_tpss_set_params)

    subroutine XC(mgga_c_bc95_set_params)(p, css, copp) bind(c)
      import
      type(c_ptr), value :: p
      real(RTYPE), value :: css, copp
    end subroutine XC(mgga_c_bc95_set_params)

    subroutine XC(mgga_c_pkzb_set_params)(p, beta, d, c0_0, c0_1, c0_2, c0_3) bind(c)
      import
      type(c_ptr), value :: p
      real(RTYPE), value :: beta, d, c0_0, c0_1, c0_2, c0_3
    end subroutine XC(mgga_c_pkzb_set_params)
  end interface

  
contains

  
  !----------------------------------------------------------------
  subroutine FXC(version)(major, minor, micro)
    integer, intent(out) :: major, minor, micro

    call XC(version)(major, minor, micro)
    
  end subroutine FXC(version)

  subroutine FXC(version_string)(version)
    character(len=*), intent(out) :: version

    type(c_ptr) :: c_version
    
    c_version = XC(version_string)()
    call c_to_f_string_ptr(c_version, version)

  end subroutine FXC(version_string)

  
  !----------------------------------------------------------------
  integer function FXC(func_info_get_number)(info) result(number)
    type(FXC(func_info_t)), intent(in) :: info

    number = XC(func_info_get_number)(info%ptr)

  end function FXC(func_info_get_number)
    
  integer function FXC(func_info_get_kind)(info) result(kind)
    type(FXC(func_info_t)), intent(in) :: info
    
    kind = XC(func_info_get_kind)(info%ptr)

  end function FXC(func_info_get_kind)

  character(len=128) function FXC(func_info_get_name)(info) result(name)
    type(FXC(func_info_t)), intent(in) :: info

    call c_to_f_string_ptr(XC(func_info_get_name)(info%ptr), name)
        
  end function FXC(func_info_get_name)

  integer function FXC(func_info_get_family)(info) result(family)
    type(FXC(func_info_t)), intent(in) :: info
    
    family = XC(func_info_get_family)(info%ptr)

  end function FXC(func_info_get_family)

  integer function FXC(func_info_get_flags)(info) result(flags)
    type(FXC(func_info_t)), intent(in) :: info

    flags = XC(func_info_get_flags)(info%ptr)

  end function FXC(func_info_get_flags)

  character(len=120) function FXC(func_info_get_refs)(info, number) result(ref)
    type(FXC(func_info_t)), intent(in)    :: info
    integer,                intent(inout) :: number ! number of the reference. Must be 0 in the first call

    type(c_ptr) :: c_ref

    c_ref = XC(func_info_get_ref)(info%ptr, number)
    
    if (c_associated(c_ref)) then
      call c_to_f_string_ptr(c_ref, ref)
      number = number + 1
    else
      number = -1
      ref = ""
    end if
    
  end function FXC(func_info_get_refs)


  !----------------------------------------------------------------
  subroutine FXC(func_init)(p, functional, nspin)
    type(FXC(func_t)),      intent(inout) :: p
    integer,                intent(in)    :: functional
    integer,                intent(in)    :: nspin

    integer :: ierr
    
    p%ptr = XC(func_alloc)()
    ierr = XC(func_init)(p%ptr, functional, nspin)

  end subroutine FXC(func_init)

  subroutine FXC(func_end)(p)
    type(FXC(func_t)), intent(inout) :: p

    call XC(func_end)(p%ptr)
    call XC(func_free)(p%ptr)

  end subroutine FXC(func_end)

  type(FXC(func_info_t)) function FXC(func_get_info)(p) result(info)
    type(FXC(func_t)), intent(in) :: p

    info%ptr = XC(func_get_info)(p%ptr)

  end function FXC(func_get_info)
  
  character(len=128) function FXC(functional_get_name)(number) result(name)     
    integer, intent(in) :: number

    call c_to_f_string_ptr(XC(functional_get_name)(number), name)

  end function FXC(functional_get_name)

  integer function FXC(functional_get_number)(func_string) result(number)
    character(len=*), intent(in) :: func_string
    
    number = XC(functional_get_number)(f_to_c_string(func_string))

  end function FXC(functional_get_number)
 
  integer function FXC(family_from_id)(id, family, number)
    integer, intent(in)                    :: id
    integer, intent(out), optional, target :: family, number

    type(c_ptr) c_family, c_number
    integer, pointer :: f_family, f_number

    if (present(family)) then
      f_family => family
      call c_f_pointer(c_family, f_family)
    else
      c_family = C_NULL_PTR
    end if
    if (present(number)) then
      f_number => number
      call c_f_pointer(c_number, f_number)
    else
      c_number = C_NULL_PTR
    end if
 
    FXC(family_from_id) = XC(family_from_id)(id, c_family, c_number)
    
  end function FXC(family_from_id)

  ! LDAs
  !----------------------------------------------------------------
  subroutine FXC(lda)(p, np, rho, zk, vrho, v2rho2, v3rho3)
    type(FXC(func_t)), intent(in)  :: p
    integer,           intent(in)  :: np
    real(RTYPE),       intent(in)  :: rho(*)
    real(RTYPE),       intent(out) :: zk(*), vrho(*), v2rho2(*), v3rho3(*)

    call XC(lda)(p%ptr, np, rho, zk, vrho, v2rho2, v3rho3)

  end subroutine FXC(lda)

  subroutine FXC(lda_exc)(p, np, rho, zk)
    type(FXC(func_t)), intent(in)  :: p
    integer,           intent(in)  :: np
    real(RTYPE),       intent(in)  :: rho(*)
    real(RTYPE),       intent(out) :: zk(*)

    call XC(lda_exc)(p%ptr, np, rho, zk)

  end subroutine FXC(lda_exc)

  subroutine FXC(lda_exc_vxc)(p, np, rho, zk, vrho)
    type(FXC(func_t)), intent(in)  :: p
    integer,           intent(in)  :: np
    real(RTYPE),       intent(in)  :: rho(*)
    real(RTYPE),       intent(out) :: zk(*), vrho(*)

    call XC(lda_exc_vxc)(p%ptr, np, rho, zk, vrho)

  end subroutine FXC(lda_exc_vxc)

  subroutine FXC(lda_vxc)(p, np, rho, vrho)
    type(FXC(func_t)), intent(in)  :: p
    integer,           intent(in)  :: np
    real(RTYPE),       intent(in)  :: rho(*)
    real(RTYPE),       intent(out) :: vrho(*)

    call XC(lda_vxc)(p%ptr, np, rho, vrho)

  end subroutine FXC(lda_vxc)

  subroutine FXC(lda_fxc)(p, np, rho, v2rho2)
    type(FXC(func_t)), intent(in)  :: p
    integer,           intent(in)  :: np
    real(RTYPE),       intent(in)  :: rho(*)
    real(RTYPE),       intent(out) :: v2rho2(*)

    call XC(lda_fxc)(p%ptr, np, rho, v2rho2)

  end subroutine FXC(lda_fxc)

  subroutine FXC(lda_kxc)(p, np, rho, v3rho3)
    type(FXC(func_t)), intent(in)  :: p
    integer,           intent(in)  :: np
    real(RTYPE),       intent(in)  :: rho(*)
    real(RTYPE),       intent(out) :: v3rho3(*)

    call XC(lda_kxc)(p%ptr, np, rho, v3rho3)

  end subroutine FXC(lda_kxc)

  subroutine FXC(lda_x_set_params)(p, alpha, relativistic, omega)
    type(FXC(func_t)), intent(inout) :: p
    real(RTYPE),       intent(in)    :: alpha, omega
    integer,           intent(in)    :: relativistic

    call XC(lda_x_set_params)(p%ptr, alpha, relativistic, omega)    
    
  end subroutine FXC(lda_x_set_params)

  subroutine FXC(lda_x_1d_set_params)(p, interaction, bb)
    type(FXC(func_t)), intent(inout) :: p
    integer,           intent(in)    :: interaction
    real(RTYPE),       intent(in)    :: bb

    call XC(lda_x_1d_set_params)(p%ptr, interaction, bb)

  end subroutine FXC(lda_x_1d_set_params)

  subroutine FXC(lda_c_1d_csc_set_params)(p, interaction, bb)
    type(FXC(func_t)), intent(inout) :: p
    integer,           intent(in)    :: interaction
    real(RTYPE),       intent(in)    :: bb

    call XC(lda_c_1d_csc_set_params)(p%ptr, interaction, bb)

  end subroutine FXC(lda_c_1d_csc_set_params)
  
  subroutine FXC(lda_c_xalpha_set_params)(p, alpha)
    type(FXC(func_t)), intent(inout) :: p
    real(RTYPE),       intent(in)    :: alpha

    call XC(lda_c_xalpha_set_params)(p%ptr, alpha)
    
  end subroutine FXC(lda_c_xalpha_set_params)
  
  subroutine FXC(lda_c_2d_prm_set_params)(p, n)
    type(FXC(func_t)), intent(inout) :: p   
    real(RTYPE),       intent(in)    :: n

    call XC(lda_c_2d_prm_set_params)(p%ptr, n)
    
  end subroutine FXC(lda_c_2d_prm_set_params)

  subroutine FXC(lda_c_vwn_set_params)(p, spin_interpolation)
    type(FXC(func_t)), intent(inout) :: p   
    integer,           intent(in)    :: spin_interpolation

    call XC(lda_c_vwn_set_params)(p%ptr, spin_interpolation)

  end subroutine FXC(lda_c_vwn_set_params)

  subroutine FXC(lda_xc_ksdt_set_params)(p, t)
    type(FXC(func_t)), intent(inout) :: p
    real(RTYPE),       intent(in)    :: t

    call XC(lda_xc_ksdt_set_params)(p%ptr, t)
    
  end subroutine FXC(lda_xc_ksdt_set_params)

  
  ! GGAs
  !----------------------------------------------------------------
  subroutine FXC(gga)(p, np, rho, sigma, zk, vrho, vsigma, &
    v2rho2, v2rhosigma, v2sigma2, v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3)
    type(FXC(func_t)), intent(in)  :: p
    integer,           intent(in)  :: np
    real(RTYPE),       intent(in)  :: rho(*), sigma(*)
    real(RTYPE),       intent(out) :: zk(*), vrho(*), vsigma(*)
    real(RTYPE),       intent(out) :: v2rho2(*), v2rhosigma(*), v2sigma2(*)
    real(RTYPE),       intent(out) :: v3rho3(*), v3rho2sigma(*), v3rhosigma2(*), v3sigma3(*)

    call XC(gga)(p%ptr, np, rho, sigma, zk, vrho, vsigma, &
      v2rho2, v2rhosigma, v2sigma2, v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3)
    
  end subroutine FXC(gga)

  subroutine FXC(gga_exc)(p, np, rho, sigma, zk)
    type(FXC(func_t)), intent(in)  :: p
    integer,           intent(in)  :: np
    real(RTYPE),       intent(in)  :: rho(*), sigma(*)
    real(RTYPE),       intent(out) :: zk(*)

    call XC(gga_exc)(p%ptr, np, rho, sigma, zk)

  end subroutine FXC(gga_exc)

  subroutine FXC(gga_exc_vxc)(p, np, rho, sigma, zk, vrho, vsigma)
    type(FXC(func_t)), intent(in)  :: p
    integer,           intent(in)  :: np
    real(RTYPE),       intent(in)  :: rho(*), sigma(*)
    real(RTYPE),       intent(out) :: zk(*), vrho(*), vsigma(*)

    call XC(gga_exc_vxc)(p%ptr, np, rho, sigma, zk, vrho, vsigma)

  end subroutine FXC(gga_exc_vxc)

  subroutine FXC(gga_vxc)(p, np, rho, sigma, vrho, vsigma)
    type(FXC(func_t)), intent(in)  :: p
    integer,           intent(in)  :: np
    real(RTYPE),       intent(in)  :: rho(*), sigma(*)
    real(RTYPE),       intent(out) :: vrho(*), vsigma(*)

    call XC(gga_vxc)(p%ptr, np, rho, sigma, vrho, vsigma)

  end subroutine FXC(gga_vxc)

  subroutine FXC(gga_fxc)(p, np, rho, sigma, v2rho2, v2rhosigma, v2sigma2)
    type(FXC(func_t)), intent(in)  :: p
    integer,           intent(in)  :: np
    real(RTYPE),       intent(in)  :: rho(*), sigma(*)
    real(RTYPE),       intent(out) :: v2rho2(*), v2rhosigma(*), v2sigma2(*)

    call XC(gga_fxc)(p%ptr, np, rho, sigma, v2rho2, v2rhosigma, v2sigma2)
    
  end subroutine FXC(gga_fxc)

  subroutine FXC(gga_kxc)(p, np, rho, sigma, v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3)
    type(FXC(func_t)), intent(in)  :: p
    integer,           intent(in)  :: np
    real(RTYPE),       intent(in)  :: rho(*), sigma(*)
    real(RTYPE),       intent(out) :: v3rho3(*), v3rho2sigma(*), v3rhosigma2(*), v3sigma3(*)

    call XC(gga_kxc)(p%ptr, np, rho, sigma, v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3)

  end subroutine FXC(gga_kxc)

  subroutine FXC(gga_lb_modified)(p, np, rho, grho, r, dedd)
    type(FXC(func_t)), intent(in)  :: p
    integer,           intent(in)  :: np
    real(RTYPE),       intent(in)  :: rho(*), grho(*)
    real(RTYPE),       intent(in)  :: r
    real(RTYPE),       intent(out) :: dedd(*)

    call XC(gga_lb_modified)(p%ptr, np, rho, grho, r, dedd)
    
  end subroutine FXC(gga_lb_modified)

  subroutine FXC(gga_x_b86_set_params)(p, beta, gamma, omega)
    type(FXC(func_t)), intent(inout) :: p
    real(RTYPE),       intent(in)    :: beta, gamma, omega

    call XC(gga_x_b86_set_params)(p%ptr, beta, gamma, omega)
    
  end subroutine FXC(gga_x_b86_set_params)

  subroutine FXC(gga_x_b88_set_params)(p, beta, gamma)
    type(FXC(func_t)), intent(inout) :: p
    real(RTYPE),       intent(in)    :: beta, gamma

    call XC(gga_x_b88_set_params)(p%ptr, beta, gamma)
    
  end subroutine FXC(gga_x_b88_set_params)

  subroutine FXC(gga_x_pbe_set_params)(p, kappa, mu)
    type(FXC(func_t)), intent(inout)  :: p
    real(RTYPE),       intent(in)     :: kappa, mu

    call XC(gga_x_pbe_set_params)(p%ptr, kappa, mu)
    
  end subroutine FXC(gga_x_pbe_set_params)

  subroutine FXC(gga_c_pbe_set_params)(p, beta)
    type(FXC(func_t)), intent(inout)  :: p
    real(RTYPE),       intent(in)     :: beta

    call XC(gga_c_pbe_set_params)(p%ptr, beta)
    
  end subroutine FXC(gga_c_pbe_set_params)

  subroutine FXC(gga_x_pw91_set_params)(p, a, b, c, d, f, alpha, expo)
    type(FXC(func_t)), intent(inout) :: p
    real(RTYPE),       intent(in)    :: a, b, c, d, f, alpha, expo

    call XC(gga_x_pw91_set_params)(p%ptr, a, b, c, d, f, alpha, expo)

  end subroutine FXC(gga_x_pw91_set_params)

  subroutine FXC(gga_x_pw91_set_params2)(p, bt, alpha, expo)
    type(FXC(func_t)), intent(inout) :: p
    real(RTYPE),       intent(in)    :: bt, alpha, expo

    call XC(gga_x_pw91_set_params2)(p%ptr, bt, alpha, expo)

  end subroutine FXC(gga_x_pw91_set_params2)

  subroutine FXC(gga_x_rpbe_set_params)(p, kappa, mu)
    type(FXC(func_t)), intent(inout) :: p
    real(RTYPE),       intent(in)    :: kappa, mu
    
    call XC(gga_x_rpbe_set_params)(p%ptr, kappa, mu)

  end subroutine FXC(gga_x_rpbe_set_params)

  subroutine FXC(gga_x_optx_set_params)(p, a, b, gamma)
    type(FXC(func_t)), intent(inout) :: p
    real(RTYPE),       intent(in)    :: a, b, gamma
    
    call XC(gga_x_optx_set_params)(p%ptr, a, b, gamma)

  end subroutine FXC(gga_x_optx_set_params)

  subroutine FXC(gga_c_lyp_set_params)(p, a, b, c, d)
    type(FXC(func_t)), intent(inout) :: p
    real(RTYPE),       intent(in)    :: a, b, c, d

    call XC(gga_c_lyp_set_params)(p%ptr, a, b, c, d)

  end subroutine FXC(gga_c_lyp_set_params)

  subroutine FXC(gga_lb_set_params)(p, modified, threshold, ip, qtot)
    type(FXC(func_t)), intent(inout) :: p
    integer,           intent(in)    :: modified
    real(RTYPE),       intent(in)    :: threshold, ip, qtot

    call XC(gga_lb_set_params)(p%ptr, modified, threshold, ip, qtot)
    
  end subroutine FXC(gga_lb_set_params)

  subroutine FXC(gga_k_tflw_set_params)(p, gamma, lambda, n)
    type(FXC(func_t)), intent(inout) :: p
    real(RTYPE),       intent(in)    :: gamma, lambda
    integer,           intent(in)    :: n

    call XC(gga_k_tflw_set_params)(p%ptr, gamma, lambda, n)

  end subroutine FXC(gga_k_tflw_set_params)
    
  subroutine FXC(gga_x_2d_b88_set_params)(p, beta)
    type(FXC(func_t)), intent(inout) :: p
    real(RTYPE),       intent(in)    :: beta

    call XC(gga_x_2d_b88_set_params)(p%ptr, beta)
    
  end subroutine FXC(gga_x_2d_b88_set_params)

  subroutine FXC(gga_x_wpbeh_set_params)(p, omega)
    type(FXC(func_t)), intent(inout) :: p
    real(RTYPE),       intent(in)    :: omega

    call XC(gga_x_wpbeh_set_params)(p%ptr, omega)

  end subroutine FXC(gga_x_wpbeh_set_params)

  subroutine FXC(gga_x_hjs_set_params)(p, omega)
    type(FXC(func_t)), intent(inout) :: p
    real(RTYPE),       intent(in)    :: omega

    call XC(gga_x_hjs_set_params)(p%ptr, omega)
    
  end subroutine FXC(gga_x_hjs_set_params)

  subroutine FXC(gga_x_ityh_set_params)(p, func_id, omega)
    type(FXC(func_t)), intent(inout) :: p
    integer,           intent(in)    :: func_id
    real(RTYPE),       intent(in)    :: omega

    call XC(gga_x_ityh_set_params)(p%ptr, func_id, omega)

  end subroutine FXC(gga_x_ityh_set_params)
  
  subroutine FXC(gga_x_sfat_set_params)(p, func_id, omega)
    type(FXC(func_t)), intent(inout) :: p
    integer,           intent(in)    :: func_id
    real(RTYPE),       intent(in)    :: omega

    call XC(gga_x_sfat_set_params)(p%ptr, func_id, omega)

  end subroutine FXC(gga_x_sfat_set_params)
    
  subroutine FXC(gga_x_ssb_sw_set_params)(p, a, b, c, d, e)
    type(FXC(func_t)), intent(inout) :: p
    real(RTYPE),       intent(in)    :: a, b, c, d, e

    call XC(gga_x_ssb_sw_set_params)(p%ptr, a, b, c, d, e)

  end subroutine FXC(gga_x_ssb_sw_set_params)

  subroutine FXC(gga_x_kt_set_params)(p, gamma, delta)
    type(FXC(func_t)), intent(inout) :: p
    real(RTYPE),       intent(in)    :: gamma, delta

    call XC(gga_x_kt_set_params)(p%ptr, gamma, delta)

  end subroutine FXC(gga_x_kt_set_params)

  subroutine FXC(gga_x_lambda_set_params)(p, n)
    type(FXC(func_t)), intent(inout) :: p
    integer,           intent(in)    :: n

    call XC(gga_x_lambda_set_params)(p%ptr, n)
    
  end subroutine FXC(gga_x_lambda_set_params)

  real(RTYPE) function FXC(gga_ak13_get_asymptotic)(homo) result(asymptotic)
    real(RTYPE), intent(in) :: homo

    asymptotic = XC(gga_ak13_get_asymptotic)(homo)
    
  end function FXC(gga_ak13_get_asymptotic)

  real(RTYPE) function FXC(hyb_exx_coef)(p) result(coef)
    type(FXC(func_t)), intent(in) :: p

    coef = XC(hyb_exx_coef)(p%ptr)
    
  end function FXC(hyb_exx_coef)
  
  subroutine FXC(hyb_cam_coef)(p, omega, alpha, beta)
    type(FXC(func_t)), intent(in)  :: p
    real(RTYPE),       intent(out) :: omega, alpha, beta

    call XC(hyb_cam_coef)(p%ptr, omega, alpha, beta)

  end subroutine FXC(hyb_cam_coef)
  
  subroutine FXC(nlc_coef)(p, nlc_b, nlc_c)
    type(FXC(func_t)), intent(in)  :: p
    real(RTYPE),       intent(out) :: nlc_b, nlc_c

    call XC(nlc_coef)(p%ptr, nlc_b, nlc_c)
    
  end subroutine FXC(nlc_coef)
  
  subroutine FXC(hyb_gga_xc_hse_set_params)(p, alpha, omega)
    type(FXC(func_t)), intent(inout) :: p
    real(RTYPE),       intent(in)    :: alpha, omega

    call XC(hyb_gga_xc_hse_set_params)(p%ptr, alpha, omega)
    
  end subroutine FXC(hyb_gga_xc_hse_set_params)
  
  subroutine FXC(hyb_gga_xc_pbeh_set_params)(p, alpha)
    type(FXC(func_t)), intent(inout) :: p
    real(RTYPE),       intent(in)    :: alpha
    
    call XC(hyb_gga_xc_pbeh_set_params)(p%ptr, alpha)

  end subroutine FXC(hyb_gga_xc_pbeh_set_params)


  ! the meta-GGAs
  !----------------------------------------------------------------
  subroutine FXC(mgga)(p, np, rho, sigma, lapl, tau, zk, vrho, vsigma, vlapl, vtau, &
    v2rho2, v2sigma2, v2lapl2, v2tau2, v2rhosigma, v2rholapl, v2rhotau, &
    v2sigmalapl, v2sigmatau, v2lapltau)
    type(FXC(func_t)), intent(in)  :: p
    integer,           intent(in)  :: np
    real(RTYPE),       intent(in)  :: rho(*), sigma(*), lapl(*), tau(*)
    real(RTYPE),       intent(out) :: zk(*), vrho(*), vsigma(*), vlapl(*), vtau(*)
    real(RTYPE),       intent(out) :: v2rho2(*), v2sigma2(*), v2lapl2(*), v2tau2(*), v2rhosigma(*), v2rholapl(*), &
                                      v2rhotau(*), v2sigmalapl(*), v2sigmatau(*), v2lapltau(*)

    call XC(mgga)(p%ptr, np, rho, sigma, lapl, tau, zk, vrho, vsigma, vlapl, vtau, &
      v2rho2, v2sigma2, v2lapl2, v2tau2, v2rhosigma, v2rholapl, v2rhotau, &
      v2sigmalapl, v2sigmatau, v2lapltau)
    
  end subroutine FXC(mgga)

  subroutine FXC(mgga_exc)(p, np, rho, sigma, lapl, tau, zk)
    type(FXC(func_t)), intent(in)  :: p
    integer,           intent(in)  :: np
    real(RTYPE),       intent(in)  :: rho(*), sigma(*), lapl(*), tau(*)
    real(RTYPE),       intent(out) :: zk(*)

    call XC(mgga_exc)(p%ptr, np, rho, sigma, lapl, tau, zk)
    
  end subroutine FXC(mgga_exc)

  subroutine FXC(mgga_exc_vxc)(p, np, rho, sigma, lapl, tau, zk, vrho, vsigma, vlapl, vtau)
    type(FXC(func_t)), intent(in)  :: p
    integer,           intent(in)  :: np
    real(RTYPE),       intent(in)  :: rho(*), sigma(*), lapl(*), tau(*)
    real(RTYPE),       intent(out) :: zk(*), vrho(*), vsigma(*), vlapl(*), vtau(*)

    call XC(mgga_exc_vxc)(p%ptr, np, rho, sigma, lapl, tau, zk, vrho, vsigma, vlapl, vtau)
    
  end subroutine FXC(mgga_exc_vxc)

  subroutine FXC(mgga_vxc)(p, np, rho, sigma, lapl, tau, vrho, vsigma, vlapl, vtau)
    type(FXC(func_t)), intent(in)  :: p
    integer,           intent(in)  :: np
    real(RTYPE),       intent(in)  :: rho(*), sigma(*), lapl(*), tau(*)
    real(RTYPE),       intent(out) :: vrho(*), vsigma(*), vlapl(*), vtau(*)

    call XC(mgga_vxc)(p%ptr, np, rho, sigma, lapl, tau, vrho, vsigma, vlapl, vtau)

  end subroutine FXC(mgga_vxc)

  subroutine FXC(mgga_fxc)(p, np, rho, sigma, lapl, tau, &
    v2rho2, v2sigma2, v2lapl2, v2tau2, v2rhosigma, v2rholapl, v2rhotau, &
    v2sigmalapl, v2sigmatau, v2lapltau)
    type(FXC(func_t)), intent(in)  :: p
    integer,           intent(in)  :: np
    real(RTYPE),       intent(in)  :: rho(*), sigma(*), lapl(*), tau(*)
    real(RTYPE),       intent(out) :: v2rho2(*), v2sigma2(*), v2lapl2(*), v2tau2(*), v2rhosigma(*), &
                                      v2rholapl(*), v2rhotau(*), v2sigmalapl(*), v2sigmatau(*), v2lapltau(*)

    call XC(mgga_fxc)(p%ptr, np, rho, sigma, lapl, tau, &
      v2rho2, v2sigma2, v2lapl2, v2tau2, v2rhosigma, v2rholapl, v2rhotau, &
      v2sigmalapl, v2sigmatau, v2lapltau)

  end subroutine FXC(mgga_fxc)

  subroutine FXC(mgga_x_tb09_set_params)(p, c)
    type(FXC(func_t)), intent(inout) :: p
    real(RTYPE),       intent(in)    :: c

    call XC(mgga_x_tb09_set_params)(p%ptr, c)

  end subroutine FXC(mgga_x_tb09_set_params)

  subroutine FXC(mgga_x_tpss_set_params)(p, b, c, e, kappa, mu)
    type(FXC(func_t)), intent(inout) :: p
    real(RTYPE),       intent(in)    :: b, c, e, kappa, mu

    call XC(mgga_x_tpss_set_params)(p%ptr, b, c, e, kappa, mu)

  end subroutine FXC(mgga_x_tpss_set_params)

  subroutine FXC(mgga_c_bc95_set_params)(p, css, copp)
    type(FXC(func_t)), intent(inout) :: p
    real(RTYPE),       intent(in)    :: css, copp

    call XC(mgga_c_bc95_set_params)(p%ptr, css, copp)

  end subroutine FXC(mgga_c_bc95_set_params)

  subroutine FXC(mgga_c_pkzb_set_params)(p, beta, d, c0_0, c0_1, c0_2, c0_3)
    type(FXC(func_t)), intent(inout) :: p
    real(RTYPE),       intent(in)    :: beta, d, c0_0, c0_1, c0_2, c0_3

    call XC(mgga_c_pkzb_set_params)(p%ptr, beta, d, c0_0, c0_1, c0_2, c0_3)

  end subroutine FXC(mgga_c_pkzb_set_params)

  
  ! Helper functions to convert between C and Fortran strings
  ! Based on the routines by Joseph M. Krahn
  function f_to_c_string(f_string) result(c_string)
    character(len=*), intent(in) :: f_string
    character(kind=c_char,len=1) :: c_string(len_trim(f_string)+1)
      
    integer :: i, strlen

    strlen = len_trim(f_string)

    forall (i=1:strlen)
      c_string(i) = f_string(i:i)
    end forall
    c_string(strlen+1) = C_NULL_CHAR

  end function f_to_c_string

  subroutine c_to_f_string(c_string, f_string)
    character(kind=c_char,len=1), intent(in)  :: c_string(*)
    character(len=*),             intent(out) :: f_string

    integer :: i

    i = 1
    do while(c_string(i) /= C_NULL_CHAR .and. i <= len(f_string))
      f_string(i:i) = c_string(i)
      i = i + 1
    end do
    if (i < len(f_string)) f_string(i:) = ' '

  end subroutine c_to_f_string

  subroutine c_to_f_string_ptr(c_string, f_string)
    type(c_ptr),      intent(in)  :: c_string
    character(len=*), intent(out) :: f_string

    character(len=1, kind=c_char), pointer :: p_chars(:)
    integer :: i

    if (.not. c_associated(c_string)) then
      f_string = ' '
    else
      call c_f_pointer(c_string, p_chars, [huge(0)])
      i = 1
      do while(p_chars(i) /= C_NULL_CHAR .and. i <= len(f_string))
        f_string(i:i) = p_chars(i)
        i = i + 1
      end do
      if (i < len(f_string)) f_string(i:) = ' '
    end if

  end subroutine c_to_f_string_ptr 
  
end module FXC(lib_m)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
