(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: work_gga_x *)

ad  := 1e-8:
b   := 0.024974:
a2  := (ad + 0.1234)/b:
a4  := 29.790:
a6  := 22.417:
a8  := 12.119:
a10 := 1570.1:
a12 := 55.944:

f0 := s-> 1 + a2*s^2 + a4*s^4 + a6*s^6 + a8*s^8 + a10*s^10 + a12*s^12:
f1 := s-> f0(s)^b/(1 + ad*s^2):

f  := x->f1(X2S*x): 