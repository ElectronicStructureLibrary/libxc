(*
 Copyright (C) 2017 M.A.L. Marques
               2022 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)
(* prefix:
  mgga_xc_gas22_params *params;

  assert(p->params != NULL);
  params = (mgga_xc_gas22_params * )(p->params);
*)

gas22_par_nx := 3:
(* More accurate value from the ipynb notebook *)
gas22_gamma_x := 0.003840616724010807:
gas22_par_x := [
    [  params_a_c_x[1], 0, 0],
    [  params_a_c_x[2], 0, 1],
    [  params_a_c_x[3], 1, 0],
]:

gas22_par_nss := 5:
(* More accurate value from the ipynb notebook *)
gas22_gamma_ss := 0.46914023462026644:
gas22_par_ss := [
    [  params_a_c_ss[1], 1, 0],
    [  params_a_c_ss[2], 0, 1],
    [  params_a_c_ss[3], 0, 2],
    [  params_a_c_ss[4], 6, 0],
    [  params_a_c_ss[5], 6, 4],
]:

gas22_par_nos := 5:
gas22_par_os := [
    [ params_a_c_os[1], 0, 0],
    [ params_a_c_os[2], 0, 2],
    [ params_a_c_os[3], 0, 6],
    [ params_a_c_os[4], 2/3, 6],
    [ params_a_c_os[5], 2/3, 2]
]:

$define lda_x_params
$include "lda_x.mpl"

$define lda_c_pw_params
$define lda_c_pw_modified_params
$include "lda_c_pw.mpl"

gas22_ux    := (mgamma, x) -> mgamma*x^2/(1 + mgamma*x^2):

(* article uses t = 2 tau convention *)
gas22_wx_ss := (t, dummy) -> (K_FACTOR_C - t)/(K_FACTOR_C + t):
gas22_wx_os := (ts0, ts1) -> (K_FACTOR_C*(ts0 + ts1) - 2*ts0*ts1)/(K_FACTOR_C*(ts0 + ts1) + 2*ts0*ts1):

gas22_gxss := (mgamma, wx, cc, n, xs, ts0, ts1) ->
  add(cc[i][1]*wx(ts0, ts1)^cc[i][2]*gas22_ux(mgamma, xs)^cc[i][3], i=1..n):
gas22_gos := (mgamma, wx, cc, n, xs, ts0, ts1) ->
  add(cc[i][1]*wx(ts0, ts1)^cc[i][2]*xs^cc[i][3], i=1..n):

gas22_fpar  := (rs, z, xs0, xs1, ts0, ts1) ->
  + lda_stoll_par(f_pw    , rs,  z,  1) * gas22_g(gas22_gamma_ss, gas22_wx_ss, gas22_par_ss, gas22_par_nss, xs0, ts0, 0)
  + lda_stoll_par(f_pw    , rs, -z, -1) * gas22_g(gas22_gamma_ss, gas22_wx_ss, gas22_par_ss, gas22_par_nss, xs1, ts1, 0):

gas22_fos := (rs, z, xs0, xs1, ts0, ts1) ->
  lda_stoll_perp(f_pw, rs, z)
  * gas22_gos(gas22_gamma_os, gas22_wx_os, gas22_par_os, gas22_par_nos, sqrt(xs0^2 + xs1^2)/sqrt(2), ts0, ts1):

gas22_f :=  (rs, z, xs0, xs1, ts0, ts1) ->
  + gas22_fpar(rs, z, xs0, xs1, ts0, ts1)
  + gas22_fos(rs, z, xs0, xs1, ts0, ts1):

gas22_f_aux := (rs, z, xs0, xs1, ts0, ts1) ->
  + opz_pow_n( z,1)/2 * f_lda_x(rs*(2/(1 + z))^(1/3),  1)
    * gas22_g(gas22_gamma_x,  gas22_wx_ss, gas22_par_x,  gas22_par_nx, xs0, ts0, 0)
  + opz_pow_n(-z,1)/2 * f_lda_x(rs*(2/(1 - z))^(1/3),  1)
    * gas22_g(gas22_gamma_x,  gas22_wx_ss, gas22_par_x,  gas22_par_nx, xs1, ts1, 0):

f :=  (rs, z, xt, xs0, xs1, us0, us1, ts0, ts1) ->
  + gas22_f_aux(rs, z, xs0, xs1, ts0, ts1)
  + gas22_f(rs, z, xs0, xs1, ts0, ts1):
