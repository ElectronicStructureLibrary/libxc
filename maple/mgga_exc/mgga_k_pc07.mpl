(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: work_mgga_x *)

a := 0.5389:
b := 3:

p := x -> X2S^2*x^2:
q := u -> X2S^2*u:

(* Equation (15) *)
fab := z -> piecewise( \
    z<=0, 0, \
    z>=a, 1, \
    (1 + exp(a/(a-z)))^b/(exp(a/z) + exp(a/(a-z)))^b \
):

(* Equation (7) *)
mDelta := (x, u) -> 8*q(u)^2/81 - p(x)*q(u)/9 + 8*p(x)^2/243:

f_W    := x -> 5*p(x)/3:

(* Equation (8) *)
f_GE4  := (x, u) -> 1 + 5*p(x)/27 + 20*q(u)/9 + mDelta(x, u):

(* Equation (11) *)
f_GE4_M := (x, u) -> f_GE4(x, u)/sqrt(1 + mDelta(x, u)^2/(1 + f_W(x))^2):

(* Equation (17) *)
f   := (rs, x, t, u) -> f_W(x) + (f_GE4_M(x, u) - f_W(x))*fab(f_GE4_M(x, u) - f_W(x)):
