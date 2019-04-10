(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)
(* replace: "dilog\(" -> "xc_dilogarithm(" *)

gg99_a := 3^(1/4)/(2*sqrt(2)*Pi^(3/2)):
gg99_b := sqrt(48)*Pi^3:

(* This is the solution of x = 2*Pi*sinh(r)/(3*cosh(r))^(1/3) *)
gg99_r := x ->
  arcsinh(sqrt(2*gg99_a^2*x^3 *cos(arccos(gg99_b/x^3)/3))):

gg99_f0 := r ->
  (-Pi^2 + 12*r*log(1 + exp(-2*r)) - 12*dilog(-exp(-2*r))) * cosh(r)^(2/3) /
  (X_FACTOR_C*2*3^(1/3)*Pi*r):

gg99_f := x -> gg99_f0(gg99_r(x)):

f := (rs, zeta, xt, xs0, xs1) ->
  gga_exchange(gg99_f, rs, zeta, xs0, xs1):
