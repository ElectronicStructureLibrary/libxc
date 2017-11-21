(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: work_gga_x *)
(* prefix:
  gga_x_sogga11_params *params;
 
  assert(p->params != NULL);
  params = (gga_x_sogga11_params * )(p->params);
*)

alpha := params_a_mu*X2S*X2S/params_a_kappa:

f0 := x -> 1 - 1/(1 + alpha*x^2):
f1 := x -> 1 - exp(-alpha*x^2):

f  := x -> add(params_a_a[i]*f0(x)^(i-1), i=1..6) + add(params_a_b[i]*f1(x)^(i-1), i=1..6):