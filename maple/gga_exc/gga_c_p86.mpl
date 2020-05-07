(*
 Copyright (C) 2017 M.A.L. Marques
               2020 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)
(* prefix:
  gga_c_p86_params *params;

  assert(p->params != NULL);
  params = (gga_c_p86_params * )(p->params);
*)

$define lda_c_pz_params
$include "lda_c_pz.mpl"

(* Equation (4) *)
DD := z  -> sqrt((1 + z)^(5/3) + (1 - z)^(5/3))/sqrt(2):

(* Equation (6) *)
CC := rs ->
  + params_a_aa
  + (params_a_bb + params_a_malpha*rs + params_a_mbeta*rs^2)/(1 + params_a_mgamma*rs + params_a_mdelta*rs^2 + 1.0e4*params_a_mbeta*rs^3):
CCinf := params_a_aa + params_a_bb:

(* Equation (9) *)
x1   := (rs, xt) -> xt/sqrt(rs/RS_FACTOR):
mPhi := (rs, xt) -> params_a_ftilde*(CCinf/CC(rs))*x1(rs, xt):

(* Equation (8) *)
H := (rs, z, xt) -> x1(rs, xt)^2*exp(-mPhi(rs, xt))*CC(rs)/DD(z):

f := (rs, z, xt, xs0, xs1) ->
  f_pz(rs, z) + H(rs, z, xt):
