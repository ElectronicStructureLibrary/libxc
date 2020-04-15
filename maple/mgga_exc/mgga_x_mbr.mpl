(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)

(* prefix:
  mgga_x_mbr_params *params;

  assert(p->params != NULL);
  params = (mgga_x_mbr_params * ) (p->params);
*)

(* replace: "br89_x\(" -> "xc_mgga_x_br89_get_x(" *)

$include "mgga_x_br89.mpl"

params_a_at := 0:

(* the three equations below are copied from mgga_x_tm.mpl *)
tm_p  := x -> (X2S*x)^2:
tm_y  := x -> (2*params_a_lambda - 1)^2 * tm_p(x):
tm_f0 := x -> (1 + 10*(70*tm_y(x)/27) + params_a_beta*tm_y(x)^2)^(1/10):

mbr_D := (ts, xs) -> 2*ts - (2*params_a_lambda - 1)^2/4 * xs^2:

br89_Q := (x, u, t) -> (
  + 6*(params_a_lambda^2 - params_a_lambda + 1/2)*(2*t - K_FACTOR_C - x^2/36)
  + 6/5*(6*Pi^2)^(2/3)*(tm_f0(x)^2 - 1)
  - 2*params_a_gamma*mbr_D(t, x)
  )/6:
