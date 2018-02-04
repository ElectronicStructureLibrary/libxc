(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: work_gga_x *)
(* prefix:
  gga_x_ssb_sw_params *params;
 
  assert(p->params != NULL);
  params = (gga_x_ssb_sw_params * )(p->params);
*)

f0 := s -> params_a_A 
   + params_a_B*s^2/(1 + params_a_C*s^2)
   - params_a_D*s^2/(1 + params_a_E*s^4):

f  := x -> f0(X2S*x):