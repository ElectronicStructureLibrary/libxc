(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: work_mgga_c *)
(* prefix:
  mgga_c_m05_params *params;

  assert(p->params != NULL);
  params = (mgga_c_m05_params * )(p->params);
*)

$define lda_c_pw_params
$define lda_c_pw_modified_params
$include "lda_c_pw.mpl"

$include "b97.mpl"

m05_comp := (rs, z, spin, xs, ts) ->
  + lda_stoll_par(f_pw, rs,  z,  1)
  * b97_g(params_a_gamma_ss, params_a_css, xs)
  * Fermi_D(xs, ts):

(* The parallel and perpendicular components of the energy *)
m05_fpar  := (rs, z, xs0, xs1, ts0, ts1) ->
  + m05_comp(rs,  z,  1, xs0, ts0)
  + m05_comp(rs, -z, -1, xs1, ts1):

m05_fperp := (rs, z, xs0, xs1, ts0, ts1) ->
  + lda_stoll_perp(f_pw, rs,  z)
  * b97_g(params_a_gamma_ab, params_a_cab, sqrt(xs0^2 + xs1^2)):

f_m05 := (rs, z, xs0, xs1, ts0, ts1) ->
  + m05_fpar (rs, z, xs0, xs1, ts0, ts1)
  + m05_fperp(rs, z, xs0, xs1, ts0, ts1):

f := (rs, z, xt, xs0, xs1, ts0, ts1, us0, us1) ->
  f_m05(rs, z, xs0, xs1, ts0, ts1):

