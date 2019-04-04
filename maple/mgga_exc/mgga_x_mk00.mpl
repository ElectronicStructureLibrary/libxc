(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)

mk00_a1 := 3*Pi/X_FACTOR_C:

mk00_f := (x, u, t) -> mk00_a1/(2*t - u/4):

f := (rs, z, xt, xs0, xs1, u0, u1, t0, t1) ->
  mgga_exchange(mk00_f, rs, z, xs0, xs1, u0, u1, t0, t1):