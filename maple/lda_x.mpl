(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: work_lda *)
(* prefix:
  lda_x_params *params;

  assert(p->params != NULL);
  params = (lda_x_params * )(p->params);
*)

$ifdef lda_x_params
params_a_alpha := 1:
$endif

lda_x_ax := -params_a_alpha*RS_FACTOR*X_FACTOR_C/2^(4/3):

f_lda_x := (rs, z) -> lda_x_ax*((1 + z)^(4/3) + (1 - z)^(4/3))/rs:
f       := (rs, z) -> f_lda_x(rs, z):