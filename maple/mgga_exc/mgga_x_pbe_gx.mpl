(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: work_mgga_x *)

$include "mgga_x_gx.mpl"

mmu := 0.001015549:
pbe_gx := x -> 1/(1 + mmu*x^2):

f := (rs, x, t, u) ->
  f_gx_a(malpha(x, t)) * pbe_gx(x):
