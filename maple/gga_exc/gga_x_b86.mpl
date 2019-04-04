(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)
(* prefix:
  gga_x_b86_params *params;

  assert(p->params != NULL);
  params = (gga_x_b86_params * )(p->params);
*)

b86_f := x -> 1 + params_a_beta*x^2/(1 + params_a_gamma*x^2)^params_a_omega:


f := (rs, zeta, xt, xs0, xs1) -> gga_exchange(b86_f, rs, zeta, xs0, xs1):
