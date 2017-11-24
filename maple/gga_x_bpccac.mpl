(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: work_gga_x *)

$define gga_x_pbe_tca_params
$include "gga_x_pbe.mpl"

$define gga_x_pw91_params
$include "gga_x_pw91.mpl"

malpha :=  1:
mbeta  := 19:

fab := x -> 1/(1 + exp(-malpha*(x - mbeta))):
f   := x -> (1 - fab(x))*f_pbe(x) + fab(x)*f_pw91(x):