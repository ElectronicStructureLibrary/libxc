(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: work_gga_c *)

$define gga_c_pbe_params
$include "gga_c_regtpss.mpl"

f1 := (rs, z, t) -> 1 + 4*A(rs, z, t)*t^2:
f2 := (rs, z, t) -> mbeta(rs, t)*(1 - f1(rs, z, t)^(-1/4))/(mgamma*A(rs, z, t)):