(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* Equation (7) from the paper *)

p     := x -> X2S^2*x^2:
z     := (x, t) -> x^2/(8*t):

alpha := (x, t) -> (t - x^2/8)/K_FACTOR_C:
qb    := (x, t) -> \
      9/20 * (alpha(x, t) - 1)/sqrt(1 + params_a_b*alpha(x, t)*(alpha(x, t) - 1)) \
      + 2*p(x)/3:

(* Equation (10) in all its glory *)
fxnum := (x, t) -> \
      + (MU_GE + params_a_c*z(x, t)^ff(z(x, t))/(1 + z(x, t)^2)^2)*p(x) \
      + 146/2025 * qb(x, t)^2 \
      - 73/405 * qb(x, t) * sqrt(1/2*(9/25*z(x, t)^2 + p(x)^2)) \
      + MU_GE^2/mkappa(x, t) * p(x)^2 \
      + 2*sqrt(params_a_e)*MU_GE*9/25*z(x, t)^2 \
      + params_a_e*params_a_mu*p(x)^3:

fxden := x -> \
      (1 + sqrt(params_a_e)*p(x))^2:

fx    := (x, t) -> fxnum(x, t)/fxden(x):
