(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: work_gga_x *)
(* prefix:
  gga_x_vmt_params *params;
 
  assert(p->params != NULL);
  params = (gga_x_vmt_params * )(p->params);
*)

f0_vmt := s -> 1 + params_a_mu*s^2*exp(-params_a_alpha*s^2)/(1 + params_a_mu*s^2):

f_vmt  := x -> f0_vmt(X2S*x):
f      := x -> f_vmt(x):  
