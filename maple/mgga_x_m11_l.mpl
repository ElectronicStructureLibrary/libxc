(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: work_mgga_c *)
(* prefix:
  mgga_x_m11_l_params *params;

  assert(p->params != NULL);
  params = (mgga_x_m11_l_params * ) (p->params);
*)

$include "mgga_x_m08.mpl"
$include "lda_x_erf.mpl"

f_spin := (rs, z, x, t) ->
  lda_x_ax*(1 + z)^(4/3)/rs * (
    + attenuation_erf(a_cnst*rs/(1 + z)^(1/3)) * m08_f(params_a_a, params_a_b, x, t)
    + (1 - attenuation_erf(a_cnst*rs/(1 + z)^(1/3))) * m08_f(params_a_c, params_a_d, x, t)
  ):


f_m11_l := (rs, z, xt, xs0, xs1, ts0, ts1) ->
  f_spin(rs, z, xs0, ts0) + f_spin(rs, -z, xs1, ts1):

f := (rs, z, xt, xs0, xs1, ts0, ts1, us0, us1) ->
  f_m11_l(rs, z, xt, xs0, xs1, ts0, ts1):
