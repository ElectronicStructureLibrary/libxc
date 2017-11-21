(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: work_mgga_x *)
(* prefix:
  mgga_x_ms_params *params;
 
  assert(pt->params != NULL);
  params = (mgga_x_ms_params * ) (pt->params);
*)

fa := a -> (1 - a^2)^3 / (1 + a^3 + params_a_b*a^6):
f0 := (p, c) -> 1 + params_a_kappa*(1 - params_a_kappa/(params_a_kappa + MU_GE*p + c)):

f := (rs, x, t, u) -> f0(X2S^2*x^2, 0) + \
  fa((t - x^2/8)/K_FACTOR_C)*(f0(X2S^2*x^2, params_a_c) - f0(X2S^2*x^2, 0)):
