(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: lda_exc *)
(* prefix:
  lda_xc_bn05_params *params;

  assert(p->params != NULL);
  params = (lda_xc_bn05_params * )(p->params);
*)

$include "lda_x_erf.mpl"
$define lda_c_pw_params
$define lda_c_pw_modified_params
$include "lda_c_pw.mpl"

p_a_hyb_omega_0_ := 1:

bn05_f := (rs, z) -> xerf_f(rs, z)
  + f_pw(rs, z)*params_a_bn05_A/(params_a_bn05_C0 + params_a_bn05_C1*rs + rs^2):

f := (rs, z) -> bn05_f(rs, z):
