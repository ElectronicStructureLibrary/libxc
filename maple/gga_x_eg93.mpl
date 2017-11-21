(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: work_gga_x *)

a1 := 1.647127:
a2 := 0.980118: 
a3 := 0.017399: 
b1 := 1.523671: 
b2 := 0.367229: 
b3 := 0.011282:

f0 := s -> (1 + a1*s^2 + a2*s^4 + a3*s^6)/(1 + b1*s^2 + b2*s^4 + b3*s^6):
f  := x -> f0(X2S*x):