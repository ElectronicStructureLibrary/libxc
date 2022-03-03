(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: hgga_exc *)

(* prefix:
*)

$define gga_x_b88_params
$include "gga_x_b88.mpl"

jse_f := (x, t) -> (1 - x^2/(8*t)):

(* eq 7 *)
mvs_f := (x, u, t, ex) -> jse_f(x, t)*b88_f(x) + (1 - jse_f(x, t))*ex:

f := (rs, z, xt, xs0, xs1, u0, u1, t0, t1, ex0, ex1) ->
  hgga_exchange(mvs_f, rs, z, xs0, xs1, u0, u1, t0, t1, ex0, ex1):