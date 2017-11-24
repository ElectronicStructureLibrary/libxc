(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: work_gga_x *)

theta0 := 1.0008:
theta1 := 0.1926:
theta2 := 1.8962:

f0 := s -> s^2/(1 + s)^2:
f  := x -> theta0 + f0(X2S*x)* (theta1 + f0(X2S*x) * theta2):
