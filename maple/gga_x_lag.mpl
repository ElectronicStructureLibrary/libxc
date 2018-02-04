(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: work_gga_x *)

a1 :=   0.041106:
a2 :=   2.626712:
a3 :=   0.092070:
a4 :=   0.657946:

f0 := s-> a1 * s^a2/(1 + a3 * s^a2)^a4:
f  := x-> f0(X2S*x):
