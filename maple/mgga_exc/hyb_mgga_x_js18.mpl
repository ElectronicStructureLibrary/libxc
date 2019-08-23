(*
 Copyright (C) 2019 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)

$include "mgga_x_tm.mpl"
$include "lda_x_erf.mpl"

attenuation_erf_f2 := a ->
  1 + 24*a^2*((20*a^2 - 64*a^4)*exp(-1/(4*a^2)) - 3 - 36*a^2 + 64*a^4 + 10*sqrt(Pi)*erf(1/(2*a))):

attenuation_erf_f3 := a ->
  1 + 8*a/7*(
      + (-8*a + 256*a^3 - 576*a^5 + 3849*a^7 - 122880*a^9)*exp(-1/(4*a^2))
      + 24*a^3*(-35 + 224*a^2 - 1440*a^4 + 5120*a^6)
      + 2*sqrt(Pi)*(-2 + 60*a^2)*erf(1/(2*a))):

js18_M := x -> (2*tm_lambda - 1)^2*tm_p(x):
js18_L := (x, t) ->
  (3*(tm_lambda^2 - tm_lambda + 1/2)*(t - K_FACTOR_C - x^2/72) - (t - K_FACTOR_C)
   + 7/18*(2*tm_lambda - 1)^2*x^2)/K_FACTOR_C:

js18_A := (rs, z, x) -> a_cnst*rs/(tm_f0(x)*(1 + z)^(1/3)):

js18_DME_SR := (rs, z, x, t) ->
  + attenuation_erf   (js18_A(rs, z, x))/tm_f0(x)^2
  + attenuation_erf_f2(js18_A(rs, z, x))*7*js18_L(x, t)/(9*tm_f0(x)^4)
  + attenuation_erf_f3(js18_A(rs, z, x))*245*js18_M(x)/(54*tm_f0(x)^4):

(* This expression (10) has \tilde A, and not A *)
js18_f_SR := (rs, z, x, t) -> tm_w(x, t)*js18_DME_SR(rs, z, x, t)
  + (1 - tm_w(x, t))*attenuation_erf(a_cnst*rs/(1 + z)^(1/3))*tm_fx_SC(x, t):

js18_f := (rs, z, x, u, t) -> -p_a_cam_beta*js18_f_SR(rs, z, x, t) + tm_f(x, u, t):

f := (rs, z, xt, xs0, xs1, u0, u1, t0, t1) -> mgga_exchange_nsp(js18_f, rs, z, xs0, xs1, u0, u1, t0, t1):
