(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: lda_exc *)

$define lda_x_params
$include "lda_x.mpl"

xrel_beta := rs -> (9*Pi/4)^(1/3)/(rs*M_C):
xrel_phi  := rs -> 1
  - 1.5*(sqrt(1 + xrel_beta(rs)^2)/xrel_beta(rs) - arcsinh(xrel_beta(rs))/xrel_beta(rs)^2)^2:

xrel_f := (rs, z) -> lda_x_f(rs, z)*xrel_phi(rs):

f := (rs, z) -> xrel_f(rs, z):
