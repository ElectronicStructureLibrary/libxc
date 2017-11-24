(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: work_mgga_x *)
(* prefix:
  mgga_x_tpss_params *params;
 
  assert(pt->params != NULL);
  params = (mgga_x_tpss_params * ) (pt->params);
*)

ff     := z -> params_a_BLOC_a + params_a_BLOC_b*z:
mkappa := (x, t) -> params_a_kappa:

$include "tpss_x.mpl"

(* Equation (5) *)

a1  := (x, t) -> mkappa(x, t)/(mkappa(x, t) + fx(x, t)):
f   := (rs, x, t, u) -> 1 + mkappa(x, t)*(1 - a1(x, t)):
