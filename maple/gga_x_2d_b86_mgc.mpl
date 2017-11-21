(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: work_gga_x *)

mbeta  := 0.003317:
mgamma := 0.008323:

f := x -> 1 + mbeta/X_FACTOR_C*x^2/(1 + mgamma*x^2)^(3/4):