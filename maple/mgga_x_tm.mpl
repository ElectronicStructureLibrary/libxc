(*
 Copyright (C) 2017 M.A.L. Marques
               2019 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: work_mgga_x *)

mlambda := 0.6866:
mbeta   := 79.873:

(* below Equation (6) *)
p  := x -> (X2S*x)^2:
y  := x -> (2*mlambda - 1)^2 * p(x):

(* Equation (7) *)
f0 := x -> (1 + 10*(70*y(x)/27) + mbeta*y(x)^2)^(1/10):

(* after Equation (9) *)
R  := (x, t) -> 1 + 595*(2*mlambda - 1)^2 * p(x)/54 \
   - (t - 3*(mlambda^2 - mlambda + 1/2)*(t - K_FACTOR_C - x^2/72))/K_FACTOR_C:

fx_DME := (x, t) -> 1/f0(x)^2 + 7*R(x, t)/(9*f0(x)^4):

malpha := (x, t) -> (t - x^2/8)/K_FACTOR_C:

(* after Equation (11) *)
qtilde := (x, t) -> 9/20*(malpha(x, t) - 1) + 2*p(x)/3:

(* Ratio tW/t; we have to make sure it's 1 at maximum *)
tratio := t -> m_min(1.0, K_FACTOR_C/t):

fx_SC := (x, t) -> (1 + 10*( \
       + (MU_GE + 50*p(x)/729)*p(x) + 146*qtilde(x, t)^2/2025 \
       - 73*qtilde(x,t)/405*(3/5*tratio(t))*(1 - tratio(t)))
       )^(1/10):

(* Equation 10 and below *)
w := t-> (tratio(t)^2 + 3*tratio(t)^3)/(1 + tratio(t)^3)^2:

f := (rs, x, t, u) -> w(t)*fx_DME(x, t) + (1 - w(t))*fx_SC(x, t):
