(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: work_mgga_c *)

$include "gga_c_gapc.mpl"

(* override definition of gap_C *)
gap_G := (rs, z, xt, par) -> RS_FACTOR^2/8 * xt^2/rs^2:

f_kcis := (rs, z, xt, xs0, xs1, ts0, ts1) ->
  + f_gap(rs, z, xt)
  - xs0^2/(8*ts0) * (1 + z)/2 * f_gap(rs,  1, xs0)
  - xs1^2/(8*ts1) * (1 - z)/2 * f_gap(rs, -1, xs1):

f  := (rs, z, xt, xs0, xs1, ts0, ts1, us0, us1) ->
  f_kcis(rs, z, xt, xs0, xs1, ts0, ts1):
