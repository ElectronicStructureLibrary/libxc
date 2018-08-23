(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: lda_exc *)
(* prefix:
  lda_k_gds08_params *params;

  assert(p->params != NULL);
  params = (lda_k_gds08_params * )(p->params);
*)

# Eq. (12)
f := (rs, zeta) -> params_a_A + params_a_B*log(n_total(rs))
  + params_a_C*log(n_total(rs))^2:
