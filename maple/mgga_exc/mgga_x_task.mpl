(*
 Copyright (C) 2019 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)

task_alpha := (x, t) -> (t/K_FACTOR_C) * m_max(1 - x^2/(8*t), 1e-10):

task_gx := x -> my_piecewise3(x > 0, 1 - exp(-4.9479*x^(-1/4)), 0):

task_a_coeff := [0.938719, -0.076371, -0.0150899]:
task_hx1 := r -> simplify(add(task_a_coeff[i+1]*ChebyshevT(i, (r - 1)/(r + 1)), i=0..2)):

task_b_coeff := [-0.628591, -2.10315, -0.5, 0.103153, 0.128591]:
task_fx  := r -> simplify(add(task_b_coeff[i+1]*ChebyshevT(i, (r - 1)/(r + 1)), i=0..4)):

task_h0x := 1.174:
task_d   := 10.0:

task_f0 := (s, a) -> task_h0x*task_gx(s^2) +
  (1.0 - task_fx(a))*(task_hx1(s^2) - task_h0x)*task_gx(s^2)^task_d:

task_f := (x, u, t) -> task_f0(X2S*x, task_alpha(x, t)):

f := (rs, z, xt, xs0, xs1, u0, u1, t0, t1) -> 
  mgga_exchange(task_f, rs, z, xs0, xs1, u0, u1, t0, t1):
