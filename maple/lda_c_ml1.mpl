(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: work_lda *)
(* prefix:
  lda_c_ml1_params *params;

  assert(p->params != NULL);
  params = (lda_c_ml1_params * )(p->params);
*)

C := 6.187335:
b := [2.763169, 1.757515, 1.741397, 0.568985, 1.572202, 1.885389]:

malpha := z -> params_a_fc*((1 + z)^params_a_q + (1 - z)^params_a_q):
mbeta  := z -> (1 + z)^(1/3)*(1 - z)^(1/3)/((1 + z)^(1/3) + (1 - z)^(1/3)):

k := (rs, z) -> C*malpha(z)*mbeta(z)*RS_FACTOR/rs:
Q := (rs, z) -> 
  - b[1]/(1 + b[2]*k(rs, z)) 
  + b[3]*log(1 + b[4]/k(rs, z))/k(rs, z)
  + b[5]/k(rs, z)
  - b[6]/k(rs, z)^2:

f := (rs, zeta) -> 1/2*(RS_FACTOR/rs)^3 * (1 - zeta^2)/4 * Q(rs, zeta):