(*
 Copyright (C) 2020 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)
(* prefix:
  mgga_k_csk_params *params;

  assert(p->params != NULL);
  params = (mgga_k_csk_params * )(p->params);
*)

csk_p := x -> X2S^2*x^2:
csk_q := u -> X2S^2*u:

(* Equation (21) *)
csk_z  := (p, q) -> 20/9*q - 40/27*p:

(* Equation (22) *)
csk_f0 := (p, q, z) ->  1 + 5*p/3
  + z*my_piecewise3(-z <= 1e-10, 1, 1 - exp(-1/m_max(1e-10, -z)^params_a_csk_a))^(1/params_a_csk_a):

csk_f := (x, u) -> 
  csk_f0(csk_p(x), csk_q(u), csk_z(csk_p(x), csk_q(u))):

f := (rs, z, xt, xs0, xs1, u0, u1, t0, t1) ->
  mgga_kinetic(csk_f, rs, z, xs0, xs1, u0, u1):
