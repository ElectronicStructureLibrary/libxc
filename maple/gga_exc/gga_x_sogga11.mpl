(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)
(* prefix:
  assert(p->params != NULL);
  const gga_x_sogga11_params * const params = (gga_x_sogga11_params * const)(p->params);
*)

sogga11_alpha := params_a_mu*X2S*X2S/params_a_kappa:

sogga11_f0 := x -> 1 - 1/(1 + sogga11_alpha*x^2):
sogga11_f1 := x -> 1 - exp(-sogga11_alpha*x^2):

sogga11_f  := x -> add(params_a_a[i]*sogga11_f0(x)^(i-1), i=1..6) + add(params_a_b[i]*sogga11_f1(x)^(i-1), i=1..6):

f := (rs, zeta, xt, xs0, xs1) -> gga_exchange(sogga11_f, rs, zeta, xs0, xs1):
