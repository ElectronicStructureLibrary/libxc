(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

op_a1 := 1.5214:
op_a2 := 0.5764:
op_b1 := 1.1284:
op_b2 := 0.3183:

op_beta := (rs, z, xs0, xs1) ->
  op_qab/b88_zab(1, op_f, rs, z, xs0, xs1):

f_op := (rs, z, xt, xs0, xs1) ->
  - (1 - z^2)*n_total(rs)/4.0
  * (op_a1*op_beta(rs, z, xs0, xs1) + op_a2)
  / (op_beta(rs, z, xs0, xs1)^4 + op_b1*op_beta(rs, z, xs0, xs1)^3 + op_b2*op_beta(rs, z, xs0, xs1)^2):

f := (rs, z, xt, xs0, xs1) ->
  f_op(rs, z, xt, xs0, xs1):
