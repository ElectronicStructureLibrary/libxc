(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: lda_exc *)
(* prefix:
  lda_xc_teter_params *params;

  assert(p->params != NULL);
  params = (lda_xc_teter_params * )(p->params);
*)

teter_f  := (rs, z) ->
  - add((params_a_teter_a[i] + f_zeta(z)*params_a_teter_ap[i])*rs^(i-1), i=1..4) /
    add((params_a_teter_b[i] + f_zeta(z)*params_a_teter_bp[i])*rs^i,     i=1..4):

f := (rs, z) -> teter_f(rs, z):