(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)
(* prefix:
  mgga_x_scan_params *params;

  assert(p->params != NULL);
  params = (mgga_x_scan_params * )(p->params);
*)

scan_p     := x -> X2S^2*x^2:
scan_alpha := (x, t) -> (t - x^2/8)/K_FACTOR_C:

scan_f_alpha := a -> my_piecewise3(
  a <= 1, exp(-params_a_c1*a/(1 - a)), -params_a_d*exp(params_a_c2/(1 - a))
  ):

scan_h1x := x -> 1 + params_a_k1*(1 - params_a_k1/(params_a_k1 + x)):

scan_b2 := sqrt(5913/405000):
scan_b1 := (511/13500)/(2*scan_b2):
scan_b3 := 1/2:
scan_b4 := MU_GE^2/params_a_k1 - 1606/18225 - scan_b1^2:
scan_y  := (x, a) -> MU_GE*scan_p(x) + scan_b4*scan_p(x)^2*exp(-scan_b4*scan_p(x)/MU_GE)
  + (scan_b1*scan_p(x) + scan_b2*(1 - a)*exp(-scan_b3*(1 - a)^2))^2:

scan_a1 := 4.9479:
scan_gx := x -> 1 - exp(-scan_a1/sqrt(X2S*x)):

scan_h0x := 1.174:
scan_f   := (x, u, t) -> (scan_h1x(scan_y(x, scan_alpha(x, t)))*(1 - scan_f_alpha(scan_alpha(x, t)))
  + scan_h0x*scan_f_alpha(scan_alpha(x, t)))*scan_gx(x):

f := (rs, z, xt, xs0, xs1, u0, u1, t0, t1) -> mgga_exchange(scan_f, rs, z, xs0, xs1, u0, u1, t0, t1):