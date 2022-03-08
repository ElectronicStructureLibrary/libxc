(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: lda_exc *)
(* prefix:
  assert(p->params != NULL);
  const lda_c_wigner_params * const params = (lda_c_wigner_params * const)(p->params);
*)

f := (rs, z) -> (1 - z^2)*params_a_a/(params_a_b + rs):
