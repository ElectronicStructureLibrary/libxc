(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: work_gga_x *)

a5  := 133.983631:
a6  :=   3.217063:
a7  := 136.707378:
a8  :=   3.223476:
a9  :=   2.675484:
a10 :=   3.473804:

f1 := s -> (1 - a5*s^a6 + a7*s^a8)/(1 + a9*s^a10):

$include "gga_x_lag.mpl"
unassign('f'):

f := x -> f0(X2S*x) + f1(X2S*x):
