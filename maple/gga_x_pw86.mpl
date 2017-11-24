(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: work_gga_x *)
(* prefix:
  gga_x_pw86_params *params;
 
  assert(p->params != NULL);
  params = (gga_x_pw86_params * )(p->params);
*)

$ifdef gga_x_rpw86_params
params_a_aa := 15*0.1234:
params_a_bb := 17.33:
params_a_cc := 0.163:
$endif

f0_pw86 := s -> (1 + params_a_aa*s^2 + params_a_bb*s^4 + params_a_cc*s^6)^(1/15):
f_pw86  := x -> f0_pw86(X2S*x):

f      := x -> f_pw86(x):