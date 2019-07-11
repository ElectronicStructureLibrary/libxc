(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)
(* prefix:
  mgga_x_rlda_params *params;

  assert(p->params != NULL);
  params = (mgga_x_rlda_params * )(p->params);
*)

(* The pre-factor should be 9 Pi/4 * 3/5, with an extra
   1/2 coming from the spin sum rule. However, in the GDME
   paper Eq. (30), they multiply if by 5/3, which I think is wrong.
   As such the pre-factor below differs by 5^2/3^2 from Eq. (30)
*)
rlda_a1 := (9/20) * 3*Pi * params_a_prefactor/(2 * X_FACTOR_C):

rlda_f := (x, u, t) -> my_piecewise3(2*t - u/4 > 1e-10, 
  rlda_a1/(2*t - u/4), 0):

f := (rs, z, xt, xs0, xs1, u0, u1, t0, t1) ->
  mgga_exchange(rlda_f, rs, z, xs0, xs1, u0, u1, t0, t1):
