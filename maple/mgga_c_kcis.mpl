(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: work_mgga_c *)

$include "gga_c_gapc.mpl"

(* The gap function gap_G is simply |nabla rho|^2/(8 rho^2) in KCIS *)
gap_G := (rs, z, xt, par) -> xt^2*n_total(rs)^(2/3)/8:

(* The polarized parameters are the same as the unpolarized ones *)
(* except that c_1, c_2, c_3 are multiplied by 0.7. 1.5, and 2.59 respectively *)
gap_par0[10] = 0.06483*((9*Pi)/4)^(2/3): (* this is approximately equal to 0.23878 *)
gap_par1 = gap_par0:
gap_par1[11] = 0.7:
gap_par1[12] = 1.5:
gap_par1[13] = 2.59:

f_kcis := (rs, z, xt, xs0, xs1, ts0, ts1) ->
  + f_gap(rs, z, xt)
  - xs0^2/(8*ts0) * (1 + z)/2 * f_gap(rs,  1, xs0)
  - xs1^2/(8*ts1) * (1 - z)/2 * f_gap(rs, -1, xs1):

f  := (rs, z, xt, xs0, xs1, ts0, ts1, us0, us1) ->
  f_kcis(rs, z, xt, xs0, xs1, ts0, ts1):
