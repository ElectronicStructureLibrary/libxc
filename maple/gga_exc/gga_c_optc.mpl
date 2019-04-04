(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)

$include "gga_c_pw91.mpl"

optc_c1 := 1.1015:
optc_c2 := 0.6625:

if evalb(Polarization = "ferr") then
    optc_f2 := (rs, z, xt, xs0, xs1) -> f_pw91(rs,  1, xs0, xs0, 0):
else
    optc_f2 := (rs, z, xt, xs0, xs1) ->
      + f_pw91(rs*(2/(1 + z))^(1/3),  1, xs0, xs0, 0)*(1 + z)/2
      + f_pw91(rs*(2/(1 - z))^(1/3), -1, xs1, 0, xs1)*(1 - z)/2:
end if:

f  := (rs, z, xt, xs0, xs1) ->
  + optc_c1*f_pw91(rs, z, xt, xs0, xs1) + (optc_c2 - optc_c1)*optc_f2(rs, z, xt, xs0, xs1):
