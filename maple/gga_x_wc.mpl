(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: work_gga_x *)

mu    := 0.2195149727645171:
c     := (146/2025)*(4/9) - (73/405)*(2/3) + (mu - 10/81):
kappa := 0.8040:

f0_aux := s -> kappa + MU_GE * s^2 + (mu - MU_GE) * s^2 * exp(-s^2) + log(1 + c * s^4):
f0_wc  := s -> 1 + kappa*(1 - kappa/f0_aux(s)):
f_wc   := x -> f0_wc(X2S*x):

f      := x -> f_wc(x):
