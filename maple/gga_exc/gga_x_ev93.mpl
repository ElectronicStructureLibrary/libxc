(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)

ev93_a1 := 1.647127:
ev93_a2 := 0.980118: 
ev93_a3 := 0.017399: 
ev93_b1 := 1.523671: 
ev93_b2 := 0.367229: 
ev93_b3 := 0.011282:

ev93_f0 := s -> (1 + ev93_a1*s^2 + ev93_a2*s^4 + ev93_a3*s^6)/(1 + ev93_b1*s^2 + ev93_b2*s^4 + ev93_b3*s^6):
ev93_f  := x -> ev93_f0(X2S*x):

f := (rs, z, xt, xs0, xs1) -> gga_exchange(ev93_f, rs, z, xs0, xs1):
