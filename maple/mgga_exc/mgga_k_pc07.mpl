(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)

$ifdef mgga_k_pcopt_params
pc07_a := 1.784720:
pc07_b := 0.258304:
$else
pc07_a := 0.5389:
pc07_b := 3:
$endif

pc07_p := x -> X2S^2*x^2:
pc07_q := u -> X2S^2*u:

(* Equation (15) *)
(* Redefined with decaying exponentials to avoid inf/inf situations *)
pc07_fab := z -> my_piecewise3(
    z<=0, 0, my_piecewise3(z>=pc07_a, 1,
    exp(-pc07_a*pc07_b/z) * (1+exp(-pc07_a/(pc07_a-z)))^pc07_b/(exp(-pc07_a/z) + exp(-pc07_a/(pc07_a-z)))^pc07_b)
):

(* Equation (7) *)
pc07_Delta := (x, u) ->
  8*pc07_q(u)^2/81 - pc07_p(x)*pc07_q(u)/9 + 8*pc07_p(x)^2/243:

pc07_f_W    := x -> 5*pc07_p(x)/3:

(* Equation (8) *)
pc07_GE4  := (x, u) ->
  1 + 5*pc07_p(x)/27 + 20*pc07_q(u)/9 + pc07_Delta(x, u):

(* Equation (11) *)
pc07_GE4_M := (x, u) ->
  pc07_GE4(x, u)/sqrt(1 + pc07_Delta(x, u)^2/(1 + pc07_f_W(x))^2):

(* Equation (17) *)
pc07_f := (x, u) ->
  pc07_f_W(x) + (pc07_GE4_M(x, u) - pc07_f_W(x))*pc07_fab(pc07_GE4_M(x, u) - pc07_f_W(x)):

f := (rs, z, xt, xs0, xs1, u0, u1, t0, t1) ->
  mgga_kinetic(pc07_f, rs, z, xs0, xs1, u0, u1):
