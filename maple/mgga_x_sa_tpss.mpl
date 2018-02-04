(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: work_mgga_x *)

a := 2.413:
b := 0.348:

params_a_b      := 0.40:
params_a_c      := 1.59096:
params_a_e      := 1.537:
params_a_mu     := 0.21951:

alpha           := (x, t) -> (t - x^2/8)/K_FACTOR_C:

(* Equation (8) *)

mkappa := (x, t) -> 2*Pi/(3*sqrt(5)) * \
          sqrt(alpha(x, t) + 1)/sqrt(a + log(alpha(x, t) + b)):
ff     := 2:

$include "tpss_x.mpl"

a1  := (x, t) -> mkappa(x, t)/(mkappa(x, t) + fx(x, t)):
f   := (rs, x, t, u) -> 1 + mkappa(x, t)*(1 - a1(x, t)):

