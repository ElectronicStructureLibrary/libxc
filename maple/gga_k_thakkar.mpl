(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: work_gga_x *)

f0 := x -> 1 + 0.0055*x^2/(1 + 0.0253*x*arcsinh(x)):
f1 := x -> -0.072*x/(1 + 2*4^(1/3)*x):

f := x -> f0(x) + f1(x):
