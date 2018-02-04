(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: work_gga_c *)

$define lda_c_pz_params
$include "lda_c_pz.mpl"

malpha := 0.023266:
mbeta  := 7.389e-6:
mgamma := 8.723:
mdelta := 0.472:
aa     := 0.001667:
bb     := 0.002568:
ftilde := 1.745*0.11:

(* Equation (4) *)
DD := z  -> sqrt((1 + z)^(5/3) + (1 - z)^(5/3))/sqrt(2):

(* Equation (6) *)
CC := rs -> 
  + aa 
  + (bb + malpha*rs + mbeta*rs^2)/(1 + mgamma*rs + mdelta*rs^2 + 1.0e4*mbeta*rs^3):
CCinf := aa + bb:

(* Equation (9) *)
x1   := (rs, xt) -> xt/sqrt(rs/RS_FACTOR):
mPhi := (rs, xt) -> ftilde*(CCinf/CC(rs))*x1(rs, xt):

(* Equation (8) *)
H := (rs, z, xt) -> x1(rs, xt)^2*exp(-mPhi(rs, xt))*CC(rs)/DD(z):

f := (rs, z, xt, xs0, xs1) ->
  f_pz(rs, z) + H(rs, z, xt):