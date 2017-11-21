(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: work_mgga_x *)

params_a_mu := MU_GE:
params_a_b  := 0.4:
params_a_c  := 2.14951:
params_a_e  := 1.987:

mgamma      := 0.000023:
ff          := 3:

mkappa := (x, t) ->
  1/(mgamma/params_a_mu^2 + mgamma/params_a_mu + 1):

$include "tpss_x.mpl"

(* Equation (8) *)

f   := (rs, x, t, u) -> 1 + fx(x, t)*exp(-mgamma*fx(x, t)/params_a_mu)/(1 + fx(x, t)) \
    + (1 - exp(-mgamma*fx(x, t)^2/params_a_mu^2))*(params_a_mu/fx(x, t) - 1):
