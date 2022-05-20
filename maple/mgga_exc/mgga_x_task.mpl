(*
 Copyright (C) 2019 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)
(* prefix:
  mgga_x_task_params *params;

  assert(p->params != NULL);
  params = (mgga_x_task_params * )(p->params);
*)

(* eq 23 *)
task_gx0 := x -> 1 - exp(-params_a_task_c*x^(-1/2)):
(* For x -> 0 task_gx -> 1. The value of x for which the exponential term gives DBL_EPSILON is *)
task_gx_threshold := (params_a_task_c / log(DBL_EPSILON))^2:
task_gx := x -> my_piecewise3(x >= task_gx_threshold, task_gx0(m_max(x, task_gx_threshold)), 1):

(* Chebyshev rational function, see doi:10.1016/0021-9991(87)90002-7 *)
ChebyshevT_rational := (n, x) -> ChebyshevT(n, (x-1)/(x+1)):

(* eq 30, lh *)
task_hx1 := s -> simplify(add(params_a_task_anu[i+1]*ChebyshevT_rational(i, s^2), i=0..2)):
(* eq 30, rh *)
task_fx  := a -> simplify(add(params_a_task_bnu[i+1]*ChebyshevT_rational(i, a), i=0..4)):

(* eq 21 *)
task_f0 := (s, a) -> params_a_task_h0x*task_gx(s^2) +
  (1.0 - task_fx(a))*(task_hx1(s^2) - params_a_task_h0x)*task_gx(s^2)^params_a_task_d:

(* assemble functional *)
task_f := (x, u, t) -> task_f0(X2S*x, mgga_alpha_s_positive(x, t)):

f := (rs, z, xt, xs0, xs1, u0, u1, t0, t1) ->
  mgga_exchange(task_f, rs, z, xs0, xs1, u0, u1, t0, t1):
