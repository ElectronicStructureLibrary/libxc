!! Copyright (C) 2016 Micael Oliveira
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
#  define RTYPE c_float
#else
#  define RTYPE c_double
#endif

module xc_f03_lib_m
  use, intrinsic :: iso_c_binding
  implicit none

  private
  public :: &
    ! version
    xc_f03_version, &
    xc_f03_version_string, &
    ! func_info
    xc_f03_func_info_t, &
    xc_f03_func_info_get_number, &
    xc_f03_func_info_get_kind, &
    xc_f03_func_info_get_name, &
    xc_f03_func_info_get_family, &
    xc_f03_func_info_get_references, &
    xc_f03_func_info_get_flags, &
    xc_f03_func_info_get_n_ext_params, &
    xc_f03_func_info_get_ext_params_description, &
    xc_f03_func_info_get_ext_params_default_value, &
    ! func_reference
    xc_f03_func_reference_t, &
    xc_f03_func_reference_get_ref, &
    xc_f03_func_reference_get_doi, &
    xc_f03_func_reference_get_bibtex, &
    ! func
    xc_f03_func_t, &
    xc_f03_func_init, &
    xc_f03_func_end, &
    xc_f03_func_get_info, &
    xc_f03_functional_get_name, &
    xc_f03_functional_get_number, &
    xc_f03_family_from_id, &
    xc_f03_number_of_functionals, &
    xc_f03_maximum_name_length, &
    xc_f03_available_functional_numbers, &
    xc_f03_available_functional_names, &
    xc_f03_func_set_dens_threshold, &
    xc_f03_func_set_ext_params, &
    ! lda
    xc_f03_lda, &
    xc_f03_lda_exc, &
    xc_f03_lda_exc_vxc, &
    xc_f03_lda_vxc, &
    xc_f03_lda_fxc, &
    xc_f03_lda_kxc, &
    ! gga
    xc_f03_gga, &
    xc_f03_gga_exc, &
    xc_f03_gga_exc_vxc, &
    xc_f03_gga_vxc, &
    xc_f03_gga_fxc, &
    xc_f03_gga_kxc, &
    xc_f03_gga_lb_modified, &
    xc_f03_gga_ak13_get_asymptotic, &
    xc_f03_hyb_exx_coef, &
    xc_f03_hyb_cam_coef, &
    xc_f03_nlc_coef, &
    ! mgga
    xc_f03_mgga, &
    xc_f03_mgga_exc, &
    xc_f03_mgga_exc_vxc, &
    xc_f03_mgga_vxc, &
    xc_f03_mgga_fxc


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
    XC_FLAGS_HAVE_EXC        =     1,   &
    XC_FLAGS_HAVE_VXC        =     2,   &
    XC_FLAGS_HAVE_FXC        =     4,   &
    XC_FLAGS_HAVE_KXC        =     8,   &
    XC_FLAGS_HAVE_LXC        =    16,   &
    XC_FLAGS_1D              =    32,   &
    XC_FLAGS_2D              =    64,   &
    XC_FLAGS_3D              =   128,   &
    XC_FLAGS_HYB_CAM         =   256,   &
    XC_FLAGS_HYB_CAMY        =   512,   &
    XC_FLAGS_VV10            =  1024,   &
    XC_FLAGS_HYB_LC          =  2048,   &
    XC_FLAGS_HYB_LCY         =  4096,   &
    XC_FLAGS_STABLE          =  8192,   &
    XC_FLAGS_DEVELOPMENT     = 16384,   &
    XC_FLAGS_NEEDS_LAPLACIAN = 32768

  integer(c_int), parameter, public :: &
    XC_TAU_EXPLICIT         =     0,   &
    XC_TAU_EXPANSION        =     1

  integer(c_int), parameter, public :: &
    XC_MAX_REFERENCES       =     5

  ! List of functionals
#include "libxc_inc.f03"
  
  ! These are old names kept for compatibility
  integer(c_int), parameter, public :: &
    XC_GGA_X_BGCP           =  38,     &
    XC_GGA_C_BGCP           =  39,     &
    XC_GGA_C_VPBE           =  83,     &
    XC_GGA_XC_LB            = 160,     &
    XC_MGGA_C_CC06          = 229,     &
    XC_GGA_K_ABSR1          = 506,     &
    XC_GGA_K_ABSR2          = 507
  
  !----------------------------------------------------------------
  interface
    subroutine xc_version(major, minor, micro) bind(c)
      import
      integer(c_int), intent(out) :: major, minor, micro
    end subroutine xc_version

    type(c_ptr) function xc_version_string() bind(c)
      import
    end function xc_version_string
  end interface

  
  !----------------------------------------------------------------
  type :: xc_f03_func_info_t
    private
    type(c_ptr) :: ptr = C_NULL_PTR
  end type xc_f03_func_info_t

  interface
    integer(c_int) function xc_func_info_get_number(info) bind(c)
      import
      type(c_ptr), value :: info
    end function xc_func_info_get_number
    
    integer(c_int) function xc_func_info_get_kind(info) bind(c)
      import
      type(c_ptr), value :: info
    end function xc_func_info_get_kind

    type(c_ptr) function xc_func_info_get_name(info) bind(c)
      import
      type(c_ptr), value :: info
    end function xc_func_info_get_name

    integer(c_int) function xc_func_info_get_family(info) bind(c)
      import
      type(c_ptr), value :: info
    end function xc_func_info_get_family

    integer(c_int) function xc_func_info_get_flags(info) bind(c)
      import
      type(c_ptr), value :: info
    end function xc_func_info_get_flags

    type(c_ptr) function xc_func_info_get_references(info, number) bind(c)
      import
      type(c_ptr),    value :: info
      integer(c_int), value :: number
    end function xc_func_info_get_references

    integer(c_int) function xc_func_info_get_n_ext_params(info) bind(c)
      import
      type(c_ptr), value :: info
    end function xc_func_info_get_n_ext_params
    
    type(c_ptr) function xc_func_info_get_ext_params_description(info, number) bind(c)
      import
      type(c_ptr),    value :: info
      integer(c_int), value :: number
    end function xc_func_info_get_ext_params_description

    real(RTYPE) function xc_func_info_get_ext_params_default_value(info, number) bind(c)
      import
      type(c_ptr),    value :: info
      integer(c_int), value :: number
    end function xc_func_info_get_ext_params_default_value
  end interface

  !----------------------------------------------------------------
  type :: xc_f03_func_reference_t
    private
    type(c_ptr) :: ptr = C_NULL_PTR
  end type xc_f03_func_reference_t

  interface
    type(c_ptr) function xc_func_reference_get_ref(reference) bind(c)
      import
      type(c_ptr), value :: reference
    end function xc_func_reference_get_ref

    type(c_ptr) function xc_func_reference_get_doi(reference) bind(c)
      import
      type(c_ptr), value :: reference
    end function xc_func_reference_get_doi

    type(c_ptr) function xc_func_reference_get_bibtex(reference) bind(c)
      import
      type(c_ptr), value :: reference
    end function xc_func_reference_get_bibtex
  end interface
    
  !----------------------------------------------------------------
  type :: xc_f03_func_t
    private
    type(c_ptr) :: ptr = C_NULL_PTR
  end type xc_f03_func_t

  interface
    type(c_ptr) function xc_func_alloc() bind(c)
      import
    end function xc_func_alloc
    
    integer(c_int) function xc_func_init(p, functional, nspin) bind(c)
      import
      type(c_ptr),    value :: p
      integer(c_int), value :: functional, nspin
    end function xc_func_init

    subroutine xc_func_end(p) bind(c)
      import
      type(c_ptr), value :: p
    end subroutine xc_func_end

    subroutine xc_func_free(p) bind(c)
      import
      type(c_ptr), value :: p
    end subroutine xc_func_free

    type(c_ptr) function xc_func_get_info(p) bind(c)
      import
      type(c_ptr), value :: p
    end function xc_func_get_info
    
    type(c_ptr) function xc_functional_get_name(number) bind(c)
      import
      integer(c_int), value :: number
    end function xc_functional_get_name

    integer(c_int) function xc_functional_get_number(func_string) bind(c)
      import
      character(kind=c_char) :: func_string(*)
    end function xc_functional_get_number
 
    integer(c_int) function xc_family_from_id(id, family, number) bind(c)
      import
      integer(c_int), value :: id
      type(c_ptr),    value :: family, number
    end function xc_family_from_id

    integer(c_int) function xc_f03_number_of_functionals() bind(c, name="xc_number_of_functionals")
      import
    end function xc_f03_number_of_functionals

    integer(c_int) function xc_f03_maximum_name_length() bind(c, name="xc_maximum_name_length")
      import
    end function xc_f03_maximum_name_length

    subroutine xc_f03_available_functional_numbers(list) bind(c, name="xc_available_functional_numbers")
      import
      integer(c_int), intent(out) :: list(*)
    end subroutine xc_f03_available_functional_numbers

    subroutine xc_available_functional_names(list) bind(c)
      import
      type(c_ptr) :: list(*)
    end subroutine xc_available_functional_names
    
    subroutine xc_func_set_dens_threshold(p, dens_threshold) bind(c)
      import
      type(c_ptr), value :: p
      real(RTYPE), value :: dens_threshold
    end subroutine xc_func_set_dens_threshold

    subroutine xc_func_set_ext_params(p, ext_params) bind(c)
      import
      type(c_ptr), value      :: p
      real(RTYPE), intent(in) :: ext_params(*)
    end subroutine xc_func_set_ext_params
  end interface


  ! LDAs
  !----------------------------------------------------------------
  interface
    subroutine xc_lda(p, np, rho, zk, vrho, v2rho2, v3rho3) bind(c)
      import
      type(c_ptr),    value       :: p
      integer(c_int), value       :: np
      real(RTYPE),    intent(in)  :: rho(*)
      real(RTYPE),    intent(out) :: zk(*), vrho(*), v2rho2(*), v3rho3(*)
    end subroutine xc_lda

    subroutine xc_lda_exc(p, np, rho, zk) bind(c)
      import
      type(c_ptr),    value       :: p
      integer(c_int), value       :: np
      real(RTYPE),    intent(in)  :: rho(*)
      real(RTYPE),    intent(out) :: zk(*)
    end subroutine xc_lda_exc

    subroutine xc_lda_exc_vxc(p, np, rho, zk, vrho) bind(c)
      import
      type(c_ptr),    value       :: p
      integer(c_int), value       :: np
      real(RTYPE),    intent(in)  :: rho(*)
      real(RTYPE),    intent(out) :: zk(*), vrho(*)
    end subroutine xc_lda_exc_vxc

    subroutine xc_lda_vxc(p, np, rho, vrho) bind(c)
      import
      type(c_ptr),    value       :: p
      integer(c_int), value       :: np
      real(RTYPE),    intent(in)  :: rho(*)
      real(RTYPE),    intent(out) :: vrho(*)
    end subroutine xc_lda_vxc

    subroutine xc_lda_fxc(p, np, rho, v2rho2) bind(c)
      import
      type(c_ptr),    value       :: p
      integer(c_int), value       :: np
      real(RTYPE),    intent(in)  :: rho(*)
      real(RTYPE),    intent(out) :: v2rho2(*)
    end subroutine xc_lda_fxc

    subroutine xc_lda_kxc(p, np, rho, v3rho3) bind(c)
      import
      type(c_ptr),    value       :: p
      integer(c_int), value       :: np
      real(RTYPE),    intent(in)  :: rho(*)
      real(RTYPE),    intent(out) :: v3rho3(*)
    end subroutine xc_lda_kxc
  end interface
  

  ! GGAs
  !----------------------------------------------------------------
  interface
    subroutine xc_gga(p, np, rho, sigma, zk, vrho, vsigma, &
      v2rho2, v2rhosigma, v2sigma2, v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3) bind(c)
      import
      type(c_ptr),    value       :: p
      integer(c_int), value       :: np
      real(RTYPE),    intent(in)  :: rho(*), sigma(*)
      real(RTYPE),    intent(out) :: zk(*), vrho(*), vsigma(*)
      real(RTYPE),    intent(out) :: v2rho2(*), v2rhosigma(*), v2sigma2(*)
      real(RTYPE),    intent(out) :: v3rho3(*), v3rho2sigma(*), v3rhosigma2(*), v3sigma3(*)
    end subroutine xc_gga

    subroutine xc_gga_exc(p, np, rho, sigma, zk) bind(c)
      import
      type(c_ptr),    value       :: p
      integer(c_int), value       :: np
      real(RTYPE),    intent(in)  :: rho(*), sigma(*)
      real(RTYPE),    intent(out) :: zk(*)
    end subroutine xc_gga_exc

    subroutine xc_gga_exc_vxc(p, np, rho, sigma, zk, vrho, vsigma) bind(c)
      import
      type(c_ptr),    value       :: p
      integer(c_int), value       :: np
      real(RTYPE),    intent(in)  :: rho(*), sigma(*)
      real(RTYPE),    intent(out) :: zk(*), vrho(*), vsigma(*)      
    end subroutine xc_gga_exc_vxc

    subroutine xc_gga_vxc(p, np, rho, sigma, vrho, vsigma) bind(c)
      import
      type(c_ptr),    value       :: p
      integer(c_int), value       :: np
      real(RTYPE),    intent(in)  :: rho(*), sigma(*)
      real(RTYPE),    intent(out) :: vrho(*), vsigma(*)
    end subroutine xc_gga_vxc

    subroutine xc_gga_fxc(p, np, rho, sigma, v2rho2, v2rhosigma, v2sigma2) bind(c)
      import
      type(c_ptr),    value       :: p
      integer(c_int), value       :: np
      real(RTYPE),    intent(in)  :: rho(*), sigma(*)
      real(RTYPE),    intent(out) :: v2rho2(*), v2rhosigma(*), v2sigma2(*)
    end subroutine xc_gga_fxc

    subroutine xc_gga_kxc(p, np, rho, sigma, v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3) bind(c)
      import
      type(c_ptr),    value       :: p
      integer(c_int), value       :: np
      real(RTYPE),    intent(in)  :: rho(*), sigma(*)
      real(RTYPE),    intent(out) :: v3rho3(*), v3rho2sigma(*), v3rhosigma2(*), v3sigma3(*)      
    end subroutine xc_gga_kxc
  end interface


  interface
    subroutine xc_gga_lb_modified(p, np, rho, grho, r, dedd) bind(c)
      import
      type(c_ptr),    value       :: p
      integer(c_int), value       :: np
      real(RTYPE),    intent(in)  :: rho(*), grho(*)
      real(RTYPE),    value       :: r
      real(RTYPE),    intent(out) :: dedd(*)
    end subroutine xc_gga_lb_modified

    real(RTYPE) function xc_gga_ak13_get_asymptotic(homo) bind(c)
      import
      real(RTYPE), value :: homo
    end function xc_gga_ak13_get_asymptotic
  end interface


  interface
    real(RTYPE) function xc_hyb_exx_coef(p) bind(c)
      import
      type(c_ptr), value :: p
    end function xc_hyb_exx_coef

    subroutine xc_hyb_cam_coef(p, omega, alpha, beta) bind(c)
      import
      type(c_ptr), value       :: p
      real(RTYPE), intent(out) :: omega, alpha, beta
    end subroutine xc_hyb_cam_coef

    subroutine xc_nlc_coef(p, nlc_b, nlc_c) bind(c)
      import
      type(c_ptr), value       :: p
      real(RTYPE), intent(out) :: nlc_b, nlc_c
    end subroutine xc_nlc_coef
  end interface


  ! the meta-GGAs
  !----------------------------------------------------------------
  interface
    subroutine xc_mgga(p, np, rho, sigma, lapl, tau, zk, vrho, vsigma, vlapl, vtau, &
      v2rho2, v2sigma2, v2lapl2, v2tau2, v2rhosigma, v2rholapl, v2rhotau, &
      v2sigmalapl, v2sigmatau, v2lapltau) bind(c)
      import
      type(c_ptr),    value       :: p
      integer(c_int), value       :: np
      real(RTYPE),    intent(in)  :: rho(*), sigma(*), lapl(*), tau(*)
      real(RTYPE),    intent(out) :: zk(*), vrho(*), vsigma(*), vlapl(*), vtau(*)
      real(RTYPE),    intent(out) :: v2rho2(*), v2sigma2(*), v2lapl2(*), v2tau2(*), v2rhosigma(*), v2rholapl(*), &
                                     v2rhotau(*), v2sigmalapl(*), v2sigmatau(*), v2lapltau(*)
    end subroutine xc_mgga

    subroutine xc_mgga_exc(p, np, rho, sigma, lapl, tau, zk) bind(c)
      import
      type(c_ptr),    value       :: p
      integer(c_int), value       :: np
      real(RTYPE),    intent(in)  :: rho(*), sigma(*), lapl(*), tau(*)
      real(RTYPE),    intent(out) :: zk(*)
    end subroutine xc_mgga_exc

    subroutine xc_mgga_exc_vxc(p, np, rho, sigma, lapl, tau, zk, vrho, vsigma, vlapl, vtau) bind(c)
      import
      type(c_ptr),    value       :: p
      integer(c_int), value       :: np
      real(RTYPE),    intent(in)  :: rho(*), sigma(*), lapl(*), tau(*)
      real(RTYPE),    intent(out) :: zk(*), vrho(*), vsigma(*), vlapl(*), vtau(*)
    end subroutine xc_mgga_exc_vxc

    subroutine xc_mgga_vxc(p, np, rho, sigma, lapl, tau, vrho, vsigma, vlapl, vtau) bind(c)
      import
      type(c_ptr),    value       :: p
      integer(c_int), value       :: np
      real(RTYPE),    intent(in)  :: rho(*), sigma(*), lapl(*), tau(*)
      real(RTYPE),    intent(out) :: vrho(*), vsigma(*), vlapl(*), vtau(*)
    end subroutine xc_mgga_vxc

    subroutine xc_mgga_fxc(p, np, rho, sigma, lapl, tau, &
      v2rho2, v2sigma2, v2lapl2, v2tau2, v2rhosigma, v2rholapl, v2rhotau, &
      v2sigmalapl, v2sigmatau, v2lapltau) bind(c)
      import
      type(c_ptr),    value       :: p
      integer(c_int), value       :: np
      real(RTYPE),    intent(in)  :: rho(*), sigma(*), lapl(*), tau(*)
      real(RTYPE),    intent(out) :: v2rho2(*), v2sigma2(*), v2lapl2(*), v2tau2(*), v2rhosigma(*), v2rholapl(*), &
                                     v2rhotau(*), v2sigmalapl(*), v2sigmatau(*), v2lapltau(*)
    end subroutine xc_mgga_fxc
  end interface

  contains
  
  !----------------------------------------------------------------
  subroutine xc_f03_version(major, minor, micro)
    integer, intent(out) :: major, minor, micro

    call xc_version(major, minor, micro)
    
  end subroutine xc_f03_version

  subroutine xc_f03_version_string(version)
    character(len=*), intent(out) :: version

    type(c_ptr) :: c_version
    
    c_version = xc_version_string()
    call c_to_f_string_ptr(c_version, version)

  end subroutine xc_f03_version_string

  
  !----------------------------------------------------------------
  integer function xc_f03_func_info_get_number(info) result(number)
    type(xc_f03_func_info_t), intent(in) :: info

    number = xc_func_info_get_number(info%ptr)

  end function xc_f03_func_info_get_number
    
  integer function xc_f03_func_info_get_kind(info) result(kind)
    type(xc_f03_func_info_t), intent(in) :: info
    
    kind = xc_func_info_get_kind(info%ptr)

  end function xc_f03_func_info_get_kind

  character(len=128) function xc_f03_func_info_get_name(info) result(name)
    type(xc_f03_func_info_t), intent(in) :: info

    call c_to_f_string_ptr(xc_func_info_get_name(info%ptr), name)
        
  end function xc_f03_func_info_get_name

  integer function xc_f03_func_info_get_family(info) result(family)
    type(xc_f03_func_info_t), intent(in) :: info
    
    family = xc_func_info_get_family(info%ptr)

  end function xc_f03_func_info_get_family

  integer function xc_f03_func_info_get_flags(info) result(flags)
    type(xc_f03_func_info_t), intent(in) :: info

    flags = xc_func_info_get_flags(info%ptr)

  end function xc_f03_func_info_get_flags

  type(xc_f03_func_reference_t) function xc_f03_func_info_get_references(info, number) result(reference)
    type(xc_f03_func_info_t), intent(in)    :: info
    integer,                intent(inout) :: number ! number of the reference. Must be 0 in the first call

    type(c_ptr) :: next_ref
    
    reference%ptr = xc_func_info_get_references(info%ptr, number)
    next_ref = xc_func_info_get_references(info%ptr, number + 1)
    if (c_associated(next_ref)) then
      number = number + 1
    else
      number = -1
    end if
    
  end function xc_f03_func_info_get_references

  integer function xc_f03_func_info_get_n_ext_params(info) result(n_ext_params)
    type(xc_f03_func_info_t), intent(in) :: info
    
    n_ext_params = xc_func_info_get_n_ext_params(info%ptr)

  end function xc_f03_func_info_get_n_ext_params

  character(len=128) function xc_f03_func_info_get_ext_params_description(info, number) result(description)
    type(xc_f03_func_info_t), intent(in) :: info
    integer,                intent(in) :: number

    call c_to_f_string_ptr(xc_func_info_get_ext_params_description(info%ptr, number), description)
        
  end function xc_f03_func_info_get_ext_params_description

  real(RTYPE) function xc_f03_func_info_get_ext_params_default_value(info, number) result(val)
    type(xc_f03_func_info_t), intent(in) :: info
    integer,                intent(in) :: number

    val = xc_func_info_get_ext_params_default_value(info%ptr, number)
        
  end function xc_f03_func_info_get_ext_params_default_value

  !----------------------------------------------------------------
  character(len=120) function xc_f03_func_reference_get_ref(reference) result(ref)
    type(xc_f03_func_reference_t), intent(in) :: reference

    call c_to_f_string_ptr(xc_func_reference_get_ref(reference%ptr), ref)

  end function xc_f03_func_reference_get_ref

  character(len=120) function xc_f03_func_reference_get_doi(reference) result(doi)
    type(xc_f03_func_reference_t), intent(in) :: reference

    call c_to_f_string_ptr(xc_func_reference_get_doi(reference%ptr), doi)

  end function xc_f03_func_reference_get_doi

  character(len=120) function xc_f03_func_reference_get_bibtex(reference) result(bibtex)
    type(xc_f03_func_reference_t), intent(in) :: reference

    call c_to_f_string_ptr(xc_func_reference_get_bibtex(reference%ptr), bibtex)

  end function xc_f03_func_reference_get_bibtex


  !----------------------------------------------------------------
  subroutine xc_f03_func_init(p, functional, nspin)
    type(xc_f03_func_t),      intent(inout) :: p
    integer,                intent(in)    :: functional
    integer,                intent(in)    :: nspin

    integer :: ierr
    
    p%ptr = xc_func_alloc()
    ierr = xc_func_init(p%ptr, functional, nspin)

  end subroutine xc_f03_func_init

  subroutine xc_f03_func_end(p)
    type(xc_f03_func_t), intent(inout) :: p

    call xc_func_end(p%ptr)
    call xc_func_free(p%ptr)

  end subroutine xc_f03_func_end

  type(xc_f03_func_info_t) function xc_f03_func_get_info(p) result(info)
    type(xc_f03_func_t), intent(in) :: p

    info%ptr = xc_func_get_info(p%ptr)

  end function xc_f03_func_get_info
  
  character(len=128) function xc_f03_functional_get_name(number) result(name)     
    integer, intent(in) :: number

    call c_to_f_string_ptr(xc_functional_get_name(number), name)

  end function xc_f03_functional_get_name

  integer function xc_f03_functional_get_number(func_string) result(number)
    character(len=*), intent(in) :: func_string
    
    number = xc_functional_get_number(f_to_c_string(func_string))

  end function xc_f03_functional_get_number
 
  integer function xc_f03_family_from_id(id, family, number)
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
 
    xc_f03_family_from_id = xc_family_from_id(id, c_family, c_number)
    
  end function xc_f03_family_from_id

  subroutine xc_f03_available_functional_names(list)
    character(len=*), intent(out) :: list(*)

    integer :: n, i, maxlen
    character(kind=c_char), allocatable, target :: names(:,:)
    type(c_ptr), allocatable :: c_list(:)
    
    n = xc_f03_number_of_functionals()
    maxlen = xc_f03_maximum_name_length()

    allocate(names(maxlen, n))
    allocate(c_list(n))
    do i = 1, n
      c_list(i) = c_loc(names(1,i))
    end do
      
    call xc_available_functional_names(c_list)
    
    do i = 1, n
      call c_to_f_string_ptr(c_list(i), list(i))
    end do

    deallocate(c_list)
    deallocate(names)
    
  end subroutine xc_f03_available_functional_names

  
  subroutine xc_f03_func_set_dens_threshold(p, dens_threshold)
    type(xc_f03_func_t), intent(in) :: p
    real(RTYPE),       intent(in) :: dens_threshold

    call xc_func_set_dens_threshold(p%ptr, dens_threshold)

  end subroutine xc_f03_func_set_dens_threshold

  subroutine xc_f03_func_set_ext_params(p, ext_params)
    type(xc_f03_func_t), intent(in) :: p
    real(RTYPE),       intent(in) :: ext_params(*)

    call xc_func_set_ext_params(p%ptr, ext_params)

  end subroutine xc_f03_func_set_ext_params
  
  ! LDAs
  !----------------------------------------------------------------
  subroutine xc_f03_lda(p, np, rho, zk, vrho, v2rho2, v3rho3)
    type(xc_f03_func_t), intent(in)  :: p
    integer,           intent(in)  :: np
    real(RTYPE),       intent(in)  :: rho(*)
    real(RTYPE),       intent(out) :: zk(*), vrho(*), v2rho2(*), v3rho3(*)

    call xc_lda(p%ptr, np, rho, zk, vrho, v2rho2, v3rho3)

  end subroutine xc_f03_lda

  subroutine xc_f03_lda_exc(p, np, rho, zk)
    type(xc_f03_func_t), intent(in)  :: p
    integer,           intent(in)  :: np
    real(RTYPE),       intent(in)  :: rho(*)
    real(RTYPE),       intent(out) :: zk(*)

    call xc_lda_exc(p%ptr, np, rho, zk)

  end subroutine xc_f03_lda_exc

  subroutine xc_f03_lda_exc_vxc(p, np, rho, zk, vrho)
    type(xc_f03_func_t), intent(in)  :: p
    integer,           intent(in)  :: np
    real(RTYPE),       intent(in)  :: rho(*)
    real(RTYPE),       intent(out) :: zk(*), vrho(*)

    call xc_lda_exc_vxc(p%ptr, np, rho, zk, vrho)

  end subroutine xc_f03_lda_exc_vxc

  subroutine xc_f03_lda_vxc(p, np, rho, vrho)
    type(xc_f03_func_t), intent(in)  :: p
    integer,           intent(in)  :: np
    real(RTYPE),       intent(in)  :: rho(*)
    real(RTYPE),       intent(out) :: vrho(*)

    call xc_lda_vxc(p%ptr, np, rho, vrho)

  end subroutine xc_f03_lda_vxc

  subroutine xc_f03_lda_fxc(p, np, rho, v2rho2)
    type(xc_f03_func_t), intent(in)  :: p
    integer,           intent(in)  :: np
    real(RTYPE),       intent(in)  :: rho(*)
    real(RTYPE),       intent(out) :: v2rho2(*)

    call xc_lda_fxc(p%ptr, np, rho, v2rho2)

  end subroutine xc_f03_lda_fxc

  subroutine xc_f03_lda_kxc(p, np, rho, v3rho3)
    type(xc_f03_func_t), intent(in)  :: p
    integer,           intent(in)  :: np
    real(RTYPE),       intent(in)  :: rho(*)
    real(RTYPE),       intent(out) :: v3rho3(*)

    call xc_lda_kxc(p%ptr, np, rho, v3rho3)

  end subroutine xc_f03_lda_kxc

  
  ! GGAs
  !----------------------------------------------------------------
  subroutine xc_f03_gga(p, np, rho, sigma, zk, vrho, vsigma, &
    v2rho2, v2rhosigma, v2sigma2, v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3)
    type(xc_f03_func_t), intent(in)  :: p
    integer,           intent(in)  :: np
    real(RTYPE),       intent(in)  :: rho(*), sigma(*)
    real(RTYPE),       intent(out) :: zk(*), vrho(*), vsigma(*)
    real(RTYPE),       intent(out) :: v2rho2(*), v2rhosigma(*), v2sigma2(*)
    real(RTYPE),       intent(out) :: v3rho3(*), v3rho2sigma(*), v3rhosigma2(*), v3sigma3(*)

    call xc_gga(p%ptr, np, rho, sigma, zk, vrho, vsigma, &
      v2rho2, v2rhosigma, v2sigma2, v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3)
    
  end subroutine xc_f03_gga

  subroutine xc_f03_gga_exc(p, np, rho, sigma, zk)
    type(xc_f03_func_t), intent(in)  :: p
    integer,           intent(in)  :: np
    real(RTYPE),       intent(in)  :: rho(*), sigma(*)
    real(RTYPE),       intent(out) :: zk(*)

    call xc_gga_exc(p%ptr, np, rho, sigma, zk)

  end subroutine xc_f03_gga_exc

  subroutine xc_f03_gga_exc_vxc(p, np, rho, sigma, zk, vrho, vsigma)
    type(xc_f03_func_t), intent(in)  :: p
    integer,           intent(in)  :: np
    real(RTYPE),       intent(in)  :: rho(*), sigma(*)
    real(RTYPE),       intent(out) :: zk(*), vrho(*), vsigma(*)

    call xc_gga_exc_vxc(p%ptr, np, rho, sigma, zk, vrho, vsigma)

  end subroutine xc_f03_gga_exc_vxc

  subroutine xc_f03_gga_vxc(p, np, rho, sigma, vrho, vsigma)
    type(xc_f03_func_t), intent(in)  :: p
    integer,           intent(in)  :: np
    real(RTYPE),       intent(in)  :: rho(*), sigma(*)
    real(RTYPE),       intent(out) :: vrho(*), vsigma(*)

    call xc_gga_vxc(p%ptr, np, rho, sigma, vrho, vsigma)

  end subroutine xc_f03_gga_vxc

  subroutine xc_f03_gga_fxc(p, np, rho, sigma, v2rho2, v2rhosigma, v2sigma2)
    type(xc_f03_func_t), intent(in)  :: p
    integer,           intent(in)  :: np
    real(RTYPE),       intent(in)  :: rho(*), sigma(*)
    real(RTYPE),       intent(out) :: v2rho2(*), v2rhosigma(*), v2sigma2(*)

    call xc_gga_fxc(p%ptr, np, rho, sigma, v2rho2, v2rhosigma, v2sigma2)
    
  end subroutine xc_f03_gga_fxc

  subroutine xc_f03_gga_kxc(p, np, rho, sigma, v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3)
    type(xc_f03_func_t), intent(in)  :: p
    integer,           intent(in)  :: np
    real(RTYPE),       intent(in)  :: rho(*), sigma(*)
    real(RTYPE),       intent(out) :: v3rho3(*), v3rho2sigma(*), v3rhosigma2(*), v3sigma3(*)

    call xc_gga_kxc(p%ptr, np, rho, sigma, v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3)

  end subroutine xc_f03_gga_kxc

  subroutine xc_f03_gga_lb_modified(p, np, rho, grho, r, dedd)
    type(xc_f03_func_t), intent(in)  :: p
    integer,           intent(in)  :: np
    real(RTYPE),       intent(in)  :: rho(*), grho(*)
    real(RTYPE),       intent(in)  :: r
    real(RTYPE),       intent(out) :: dedd(*)

    call xc_gga_lb_modified(p%ptr, np, rho, grho, r, dedd)
    
  end subroutine xc_f03_gga_lb_modified

  real(RTYPE) function xc_f03_gga_ak13_get_asymptotic(homo) result(asymptotic)
    real(RTYPE), intent(in) :: homo

    asymptotic = xc_gga_ak13_get_asymptotic(homo)
    
  end function xc_f03_gga_ak13_get_asymptotic

  real(RTYPE) function xc_f03_hyb_exx_coef(p) result(coef)
    type(xc_f03_func_t), intent(in) :: p

    coef = xc_hyb_exx_coef(p%ptr)
    
  end function xc_f03_hyb_exx_coef
  
  subroutine xc_f03_hyb_cam_coef(p, omega, alpha, beta)
    type(xc_f03_func_t), intent(in)  :: p
    real(RTYPE),       intent(out) :: omega, alpha, beta

    call xc_hyb_cam_coef(p%ptr, omega, alpha, beta)

  end subroutine xc_f03_hyb_cam_coef
  
  subroutine xc_f03_nlc_coef(p, nlc_b, nlc_c)
    type(xc_f03_func_t), intent(in)  :: p
    real(RTYPE),       intent(out) :: nlc_b, nlc_c

    call xc_nlc_coef(p%ptr, nlc_b, nlc_c)
    
  end subroutine xc_f03_nlc_coef
  

  ! the meta-GGAs
  !----------------------------------------------------------------
  subroutine xc_f03_mgga(p, np, rho, sigma, lapl, tau, zk, vrho, vsigma, vlapl, vtau, &
    v2rho2, v2sigma2, v2lapl2, v2tau2, v2rhosigma, v2rholapl, v2rhotau, &
    v2sigmalapl, v2sigmatau, v2lapltau)
    type(xc_f03_func_t), intent(in)  :: p
    integer,           intent(in)  :: np
    real(RTYPE),       intent(in)  :: rho(*), sigma(*), lapl(*), tau(*)
    real(RTYPE),       intent(out) :: zk(*), vrho(*), vsigma(*), vlapl(*), vtau(*)
    real(RTYPE),       intent(out) :: v2rho2(*), v2sigma2(*), v2lapl2(*), v2tau2(*), v2rhosigma(*), v2rholapl(*), &
                                      v2rhotau(*), v2sigmalapl(*), v2sigmatau(*), v2lapltau(*)

    call xc_mgga(p%ptr, np, rho, sigma, lapl, tau, zk, vrho, vsigma, vlapl, vtau, &
      v2rho2, v2sigma2, v2lapl2, v2tau2, v2rhosigma, v2rholapl, v2rhotau, &
      v2sigmalapl, v2sigmatau, v2lapltau)
    
  end subroutine xc_f03_mgga

  subroutine xc_f03_mgga_exc(p, np, rho, sigma, lapl, tau, zk)
    type(xc_f03_func_t), intent(in)  :: p
    integer,           intent(in)  :: np
    real(RTYPE),       intent(in)  :: rho(*), sigma(*), lapl(*), tau(*)
    real(RTYPE),       intent(out) :: zk(*)

    call xc_mgga_exc(p%ptr, np, rho, sigma, lapl, tau, zk)
    
  end subroutine xc_f03_mgga_exc

  subroutine xc_f03_mgga_exc_vxc(p, np, rho, sigma, lapl, tau, zk, vrho, vsigma, vlapl, vtau)
    type(xc_f03_func_t), intent(in)  :: p
    integer,           intent(in)  :: np
    real(RTYPE),       intent(in)  :: rho(*), sigma(*), lapl(*), tau(*)
    real(RTYPE),       intent(out) :: zk(*), vrho(*), vsigma(*), vlapl(*), vtau(*)

    call xc_mgga_exc_vxc(p%ptr, np, rho, sigma, lapl, tau, zk, vrho, vsigma, vlapl, vtau)
    
  end subroutine xc_f03_mgga_exc_vxc

  subroutine xc_f03_mgga_vxc(p, np, rho, sigma, lapl, tau, vrho, vsigma, vlapl, vtau)
    type(xc_f03_func_t), intent(in)  :: p
    integer,           intent(in)  :: np
    real(RTYPE),       intent(in)  :: rho(*), sigma(*), lapl(*), tau(*)
    real(RTYPE),       intent(out) :: vrho(*), vsigma(*), vlapl(*), vtau(*)

    call xc_mgga_vxc(p%ptr, np, rho, sigma, lapl, tau, vrho, vsigma, vlapl, vtau)

  end subroutine xc_f03_mgga_vxc

  subroutine xc_f03_mgga_fxc(p, np, rho, sigma, lapl, tau, &
    v2rho2, v2sigma2, v2lapl2, v2tau2, v2rhosigma, v2rholapl, v2rhotau, &
    v2sigmalapl, v2sigmatau, v2lapltau)
    type(xc_f03_func_t), intent(in)  :: p
    integer,           intent(in)  :: np
    real(RTYPE),       intent(in)  :: rho(*), sigma(*), lapl(*), tau(*)
    real(RTYPE),       intent(out) :: v2rho2(*), v2sigma2(*), v2lapl2(*), v2tau2(*), v2rhosigma(*), &
                                      v2rholapl(*), v2rhotau(*), v2sigmalapl(*), v2sigmatau(*), v2lapltau(*)

    call xc_mgga_fxc(p%ptr, np, rho, sigma, lapl, tau, &
      v2rho2, v2sigma2, v2lapl2, v2tau2, v2rhosigma, v2rholapl, v2rhotau, &
      v2sigmalapl, v2sigmatau, v2lapltau)

  end subroutine xc_f03_mgga_fxc


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
  
end module xc_f03_lib_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
