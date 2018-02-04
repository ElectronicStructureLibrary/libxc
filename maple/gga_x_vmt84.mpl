(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: work_gga_x *)
(* prefix:
  gga_x_vmt84_params *params;
 
  assert(p->params != NULL);
  params = (gga_x_vmt84_params * )(p->params);
*)

$include "gga_x_vmt.mpl"

f0_vmt84 := s -> (1 - exp(-params_a_alpha*s^4))/s^2 - 1 + exp(-params_a_alpha*s^4):

f  := x -> f_vmt(x) + f0_vmt84(X2S*x):