(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: work_mgga_x *)

(* prefix:
  mgga_x_mvsb_params *params;

  assert(pt->params != NULL);
  params = (mgga_x_mvsb_params * ) (pt->params);
*)

$include "mgga_x_mvs.mpl"

beta := (t,x) -> alpha(t,x)*K_FACTOR_C/(t-K_FACTOR_C):

f := (rs, x, t, u) -> (1 + params_a_k0*fa(beta(t,x))) / (1 + params_a_b*(X2S*x)^4)^(1/8):
