(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

$define lda_c_pw_params
$define lda_c_pw_modified_params
$include "lda_c_pw.mpl"

ux    := (mgamma, x) -> mgamma*x^2/(1 + mgamma*x^2):

(* article uses t = 2 tau convention *)
wx_ss := (t, dummy) -> (K_FACTOR_C - t)/(K_FACTOR_C + t):
wx_os := (ts0, ts1) -> (K_FACTOR_C*(ts0 + ts1) - 2*ts0*ts1)/(K_FACTOR_C*(ts0 + ts1) + 2*ts0*ts1):

b97mv_g := (mgamma, wx, cc, n, xs, ts0, ts1) ->
  add(cc[i][1]*wx(ts0, ts1)^cc[i][2]*ux(mgamma, xs)^cc[i][3], i=1..n):

b97mv_fpar  := (lda_func, rs, z, xs0, xs1, ts0, ts1) ->
  + lda_stoll_par(lda_func, rs,  z,  1) * b97mv_g(gamma_x,  wx_ss, par_x,  par_n, xs0, ts0, 0)
  + lda_stoll_par(f_pw    , rs,  z,  1) * b97mv_g(gamma_ss, wx_ss, par_ss, par_n, xs0, ts0, 0)
  + lda_stoll_par(lda_func, rs, -z, -1) * b97mv_g(gamma_x,  wx_ss, par_x,  par_n, xs1, ts1, 0)
  + lda_stoll_par(f_pw    , rs, -z, -1) * b97mv_g(gamma_ss, wx_ss, par_ss, par_n, xs1, ts1, 0):

b97mv_fos := (rs, z, xs0, xs1, ts0, ts1) ->
  lda_stoll_perp(f_pw, rs, z)
  * b97mv_g(gamma_os, wx_os, par_os, par_n, sqrt(xs0^2 + xs1^2)/sqrt(2), ts0, ts1):

f_b97mv :=  (lda_func, rs, z, xt, xs0, xs1, ts0, ts1) ->
  + b97mv_fpar(lda_func, rs, z, xs0, xs1, ts0, ts1)
  + b97mv_fos(rs, z, xs0, xs1, ts0, ts1):
