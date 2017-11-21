(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: work_gga_c *)

$define lda_c_vbh_params
$include "lda_c_hl.mpl"

malpha := (4*Pi/3)^(1/6):

(* 4.28e-3/2, where the 2 comes from the covertion from Ryd. to Hartree *)
a1     := Pi/(16*(3*Pi^2)^(4/3)):
a2     := -0.262:
a3     := -7/(9*2^(5/3)):

t1 := (rs, z, xt) -> 
  + xt^2*exp(a2*xt/(malpha*sqrt(rs))) 
  * sqrt(2)/sqrt((1 + z)^(5/3) + (1 - z)^(5/3)):
t2 := (z, xs0, xs1) -> 
  + a3*(xs0^2*(1 + z)^(4/3) + xs1^2*(1 - z)^(4/3)):

f := (rs, z, xt, xs0, xs1) -> 
  + f_hl(rs, z)
  + a1/(malpha^2*rs)*(t1(rs, z, xt) + t2(z, xs0, xs1)):

