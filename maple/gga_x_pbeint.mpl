(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: work_gga_x *)
(* prefix:
  gga_x_pbeint_params *params;
 
  assert(p->params != NULL);
  params = (gga_x_pbeint_params * )(p->params);
*)

mu := s -> params_a_muGE + (params_a_muPBE - params_a_muGE)* \
   params_a_alpha*s^2/(1 + params_a_alpha * s^2):

(* this is the gga_x_pbe expression *)
f0 := s -> 1 + params_a_kappa * (1 - params_a_kappa/(params_a_kappa + mu(s)*s^2)):
f  := x -> f0(X2S * x):
