(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: work_lda *)
(* prefix:
  lda_xc_1d_ehwlrg_params *params;
 
  assert(p->params != NULL);
  params = (lda_xc_1d_ehwlrg_params * )(p->params);
*)

n := rs -> 1/(2*rs):

f := (rs, zeta) -> \
 (params_a_a1 + params_a_a2*n(rs) + params_a_a3*n(rs)^2) * n(rs)^params_a_alpha:
