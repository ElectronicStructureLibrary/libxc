(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: work_mgga_x *)
(* prefix:
  mgga_x_m08_params *params;

  assert(pt->params != NULL);
  params = (mgga_x_m08_params * ) (pt->params);
*)

params_a_rpbe_kappa := 0.552:
params_a_rpbe_mu    := MU_GE:
$include "gga_x_rpbe.mpl"

params_a_kappa := KAPPA_PBE:
params_a_mu    := 0.21951:
$include "gga_x_pbe.mpl"

m08_f := (a, b, x, t) ->
  + f_pbe(x) *mgga_series_w(a, 12, t)
  + f_rpbe(x)*mgga_series_w(b, 12, t):

f_m08 := (x, t) ->
  m08_f(params_a_a, params_a_b, x, t):

f := (rs, x, t, u) ->
  f_m08(x, t):
