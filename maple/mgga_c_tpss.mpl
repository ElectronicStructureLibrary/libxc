(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: work_mgga_c *)
(* prefix:
  mgga_c_tpss_params *params;

  assert(p->params != NULL);
  params = (mgga_c_tpss_params * )(p->params);
*)

(* beta is taken from the params *)
params_a_gamma := (1 - log(2))/Pi^2:
params_a_BB    := 1:
$include "gga_c_pbe.mpl"

$include "tpss.mpl"

f := (rs, z, xt, xs0, xs1, ts0, ts1, us0, us1) ->
  + f_tpss(f_pbe, rs, z, xt, xs0, xs1, ts0, ts1):




