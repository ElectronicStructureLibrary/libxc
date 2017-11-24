(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: work_gga_x *)

mbeta  :=  0.0042:
mgamma :=  6:
c0    :=  1.09878:
c1    := -2.51173:
c2    :=  0.0156233:

f_aux := x -> 1 + mgamma*mbeta*x*arcsinh(x):

f := x -> c0 + mbeta/X_FACTOR_C*x^2*(c1/f_aux(x) + c2/(mbeta*f_aux(x)^2)):