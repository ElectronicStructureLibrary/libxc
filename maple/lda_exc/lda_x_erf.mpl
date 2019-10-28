(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: lda_exc *)

$define lda_x_params
$include "lda_x.mpl"
$include "attenuation.mpl"

a_cnst := (4/(9*Pi))^(1/3)*p_a_cam_omega/2:

lda_x_erf_spin := (rs, z) -> 
  lda_x_ax*(1 + z)^(4/3)/rs * attenuation_erf(a_cnst*rs/(1 + z)^(1/3)):

if evalb(Polarization = "ferr") then
  f_lda_x_erf := (rs, z) -> lda_x_erf_spin(rs, 1)
else
  f_lda_x_erf := (rs, z) -> lda_x_erf_spin(rs, z) + lda_x_erf_spin(rs, -z):
end if:

f := (rs, z) -> f_lda_x_erf(rs, z):
