(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: work_mgga_x *)

mlambda := 0.6866:
mbeta   := 79.873:

(* below Equation (6) *)
y  := p -> (2*mlambda - 1)^2 * p:

(* Equation (7) *)
f0 := x -> (1 + 10*70*y(X2S^2*x^2)/27 + mbeta*y(X2S^2*x^2)^2)^(1/10):

R  := (x, t) -> 1 + 595*(2*mlambda - 1)^2 * X2S^2*x^2/54 \
   - (t - 3*(mlambda^2 - mlambda + 1/2)*(t - K_FACTOR_C - x^2/72))/K_FACTOR_C:

fx_DME := (x, t) -> 1/f0(x)^2 + 7*R(x, t)/(9*f0(x)^4):

malpha := (x, t) -> (t - x^2/8)/K_FACTOR_C:
qtilde := (x, t) -> 9/20*(malpha(x, t) - 1) + 2*X2S^2*x^2/3:

fx_SC_aux := (x, t) -> (1 + 10*( \
       + (MU_GE + 50*X2S^2*x^2/729)*X2S^2*x^2 \
       + 146*qtilde(x, t)^2/2025 \
       - 73*qtilde(x,t)/405*(3*K_FACTOR_C/(5*t))*(1 - K_FACTOR_C/t))
       ):

(* It turns out that the quantity fx_SC_aux can become negative, leading to NaN.
   this appears, for example, for the following densities

   ./xc-get_data 540 1 729.292719352103 0 365940004.708665 0 0 -970050.160877395 0 64600.7149099035 0

   As a turn around, I demand that the function is always larger than 1e-16
*)
fx_SC  := (x, t) -> m_max(fx_SC_aux(x,t), 1e-16)^(0.1):

w := t-> (K_FACTOR_C^2/t^2 + 3*K_FACTOR_C^3/t^3)/(1 + K_FACTOR_C^3/t^3)^2:

f := (rs, x, t, u) -> w(t)*fx_DME(x, t) + (1 - w(t))*fx_SC(x, t):
