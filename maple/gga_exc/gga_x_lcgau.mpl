(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)
(* prefix:
  gga_x_lcgau_params *params;

  assert(p->params != NULL);
  params = (gga_x_lcgau_params * )(p->params);
*)

$include "attenuation.mpl"
$define gga_x_b88_params
$include "gga_x_b88.mpl"

(* Eq. 5 - Note that their K_s = 2*X_FACTOR_C*b88_f
   Note that there is a misprint sqrt(mu) -> sqrt(pi)
*)
lcgau_arg := (rs, z, xs) ->
  sqrt(2*X_FACTOR_C*b88_f(xs)/Pi)/(6*n_spin(rs, z)^(1/3)):

(* Eq. 4 + Eq. 8 *)
lcgau_f := (rs, z, xs) ->  b88_f(xs) * (
  + attenuation_erf(params_a_omega0*lcgau_arg(rs, z, xs))
  +  params_a_alpha1*attenuation_gau(params_a_omega1*lcgau_arg(rs, z, xs))
  +  params_a_alpha2*attenuation_gau(params_a_omega2*lcgau_arg(rs, z, xs))
  ):

f := (rs, zeta, xt, xs0, xs1) -> gga_exchange_nsp(lcgau_f, rs, zeta, xs0, xs1):
