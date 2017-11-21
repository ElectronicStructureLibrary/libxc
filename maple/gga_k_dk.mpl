(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: work_gga_x *)
(* prefix:
  gga_k_dk_params *params;
 
  assert(p->params != NULL);
  params = (gga_k_dk_params * )(p->params);
*)

f := x -> 
  add(1*params_a_aa[i]*x^(2*(i-1)), i=1..5) /
  add(1*params_a_bb[i]*x^(2*(i-1)), i=1..5):