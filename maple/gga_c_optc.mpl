(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: work_gga_c *)

$include "gga_c_pw91.mpl"

c1 := 1.1015:
c2 := 0.6625:

f  := (rs, z, xt, xs0, xs1) ->
  + c1*f_pw91(rs, z, xt, xs0, xs1) + (c2 - c1)*(
  +    f_pw91(rs*(2*(1 + z))^(1/3),  1, xs0, xs0, 0)
  +    f_pw91(rs*(2*(1 - z))^(1/3), -1, xs1, 0, xs1)
):