(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: lda_exc *)
(* prefix:
  lda_xc_zlp_params *params;

  assert(p->params != NULL);
  params = (lda_xc_zlp_params * )(p->params);
*)

zlp_a0 := params_a_zlp_a0*RS_FACTOR:
zlp_k  := params_a_zlp_k *RS_FACTOR:

zlp_f := (rs, z) -> -zlp_a0*(1 - zlp_k*log(1 + rs/zlp_k)/rs)/rs:

f := (rs, z) -> zlp_f(rs, z):
