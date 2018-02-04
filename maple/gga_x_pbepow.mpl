(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: work_gga_x *)

kappa := 0.8040:
mu    := 0.2195149727645171:
m     := 100:

gamm  := m*mu/kappa:
Cx    := kappa/m:

f0 := s -> 1 + add(Cx * (gamm*s^2/(1 + gamm*s^2))^i, i=1..m):

f  := x -> f0(X2S*x):
