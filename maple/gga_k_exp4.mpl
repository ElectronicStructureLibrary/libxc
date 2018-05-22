(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: work_gga_x *)

a1 := 199.81:
a2 := 4.3476:
c1 := 0.8524:
c2 := 1.2264:

# This is Eq. (40) of the paper.
f0 := s -> c1*(1 - exp(-a1*s^2)) + c2*(1 - exp(-a2*s^4)):

f  := x -> f0(X2S*x):