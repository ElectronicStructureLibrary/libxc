(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: work_gga_x *)

a  := 0.157:
c1 := 0.21951:
c2 := -0.015:

f0 := s -> s^2/(1 + a*s^2):
f := x -> 1 + c1*f0(X2S*x) + c2*f0(X2S*x)^2: