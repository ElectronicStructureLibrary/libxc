(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: work_gga_x *)
(* prefix:
  gga_x_mpbe_params *params;
 
  assert(p->params != NULL);
  params = (gga_x_mpbe_params * )(p->params);
*)

a  := 0.157:
c1 := 0.21951:
c2 := -0.015:

f0 := s -> s^2/(1 + params_a_a*s^2):
f := x -> 1
  + params_a_c1*f0(X2S*x)
  + params_a_c2*f0(X2S*x)^2
  + params_a_c3*f0(X2S*x)^3: