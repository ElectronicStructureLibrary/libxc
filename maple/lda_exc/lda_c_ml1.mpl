(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: lda_exc *)
(* prefix:
  lda_c_ml1_params *params;

  assert(p->params != NULL);
  params = (lda_c_ml1_params * )(p->params);
*)

CC := 6.187335:
bb := [2.763169, 1.757515, 1.741397, 0.568985, 1.572202, 1.885389]:

malpha := z -> params_a_fc*(opz_pow_n(z,params_a_q) + opz_pow_n(-z,params_a_q)):
mbeta  := z -> (1 - z^2)^(1/3)/(opz_pow_n(z,1/3) + opz_pow_n(-z,1/3)):

kk := (rs, z) -> CC*malpha(z)*mbeta(z)*RS_FACTOR/rs:
QQ := (rs, z) ->
  - bb[1]/(1 + bb[2]*kk(rs, z))
  + bb[3]*log(1 + bb[4]/kk(rs, z))/kk(rs, z)
  + bb[5]/kk(rs, z)
  - bb[6]/kk(rs, z)^2:

f := (rs, zeta) -> 1/2*(RS_FACTOR/rs)^3 * (1 - zeta^2)/4 * QQ(rs, zeta):
