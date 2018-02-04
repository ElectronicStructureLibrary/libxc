(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: work_gga_x *)

Ax       := -3/4*(3/Pi)^(1/3):
mu       := 0.2195149727645171:
alpha    := -Ax*mu:
c        := alpha/(3*Pi^2)^(1/3):
alphaoAx := -mu:

f0 := s -> 1 - alphaoAx*s*log(1 + s)/(1 + c*log(1 + s)):
f  := x -> f0(X2S*x):