(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)

$define lda_c_pw_params
$include "lda_c_pw.mpl"

$define lda_x_params
$include "lda_x.mpl"

b98_gamma_x  := 0.11:
b98_gamma_ss := 1.6:
b98_gamma_ab := 0.14:

params_a_c_x  := [0.8085,  0.6682, 0.1420]:
params_a_c_ss := [0.2606, -0.9608, 0.9023]:
params_a_c_ab := [1.2033, -2.2717, 0.9596]:

b98_q := (x, u, t) -> (u/4 - t + x^2/8 + K_FACTOR_C)/K_FACTOR_C:

b98_g := (gamma, cc, q) -> add(cc[i]*(gamma*q^2/sqrt(1 + gamma*q^2))^(i-1), i=1..3):

if evalb(Polarization = "ferr") then
  b98_f := (lda_func, g_ss, cc_ss, g_ab, cc_ab, rs, z, xs0, xs1, us0, us1, ts0, ts1) ->
    + lda_func(rs,  1)
      * b98_g(g_ss, cc_ss, b98_q(xs0, us0, ts0)) * Fermi_D(xs0, ts0):
else
  b98_f := (lda_func, g_ss, cc_ss, g_ab, cc_ab, rs, z, xs0, xs1, us0, us1, ts0, ts1) ->
    + lda_stoll_par(lda_func, rs,  z,  1)
      * b98_g(g_ss, cc_ss, b98_q(xs0, us0, ts0)) * Fermi_D(xs0, ts0)
    + lda_stoll_par(lda_func, rs, -z, -1)
      * b98_g(g_ss, cc_ss, b98_q(xs1, us1, ts1)) * Fermi_D(xs1, ts1)
    + lda_stoll_perp(lda_func, rs, z)
      * b98_g(g_ab, cc_ab, (b98_q(xs0, us0, ts0) + b98_q(xs1, us1, ts1))/2):
end if:

f := (rs, z, xt, xs0, xs1, us0, us1, ts0, ts1) ->
  + f_lda_x(rs,  1)*b98_g(b98_gamma_x, params_a_c_x, b98_q(xs0, us0, ts0))
  + f_lda_x(rs, -1)*b98_g(b98_gamma_x, params_a_c_x, b98_q(xs1, us1, ts1))
  + b98_f(f_pw, b98_gamma_ss, params_a_c_ss, b98_gamma_ab, params_a_c_ab,
        rs, z, xs0, xs1, us0, us1, ts0, ts1):