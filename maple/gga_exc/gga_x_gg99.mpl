(*
 Copyright (C) 2017 M.A.L. Marques
               2019 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)
(* replace: "dilog\(" -> "xc_dilogarithm(" *)

gg99_a := 3^(1/4)/(2*sqrt(2)*Pi^(3/2)):
gg99_b := sqrt(48)*Pi^3:

(* This is the solution of x = 2*Pi*sinh(r)/(3*cosh(r))^(1/3).
*)

(* Equation 22 in the paper is

   gg99_r_b1 := x-> arcsinh( (gg99_a * sqrt(x^2 + (gg99_b + sqrt(gg99_b^2
- x^6))^(2/3))) / (gg99_b + sqrt(gg99_b^2 - x^6))^(1/6) ):

*)

gg99_r_b1 := x ->
arcsinh(sqrt(2*gg99_a^2*x^3*cos(arccos(gg99_b/x^3)/3))):

(* The second branch is from Andrew Gilbert via email *)

gg99_r_b2 := x-> arcsinh(sqrt(x^3*sqrt(3) / (4*Pi^3) *
cos(arctan(sqrt(x^6/(48*Pi^6)-1))/3))):

(* Glue the pieces together *)
gg99_r := x -> my_piecewise3( x < Pi*48^(1/6), gg99_r_b1(x), gg99_r_b2(x)):

gg99_f0 := r ->
  (-Pi^2 + 12*r*log(1 + exp(-2*r)) - 12*dilog(-exp(-2*r))) * cosh(r)^(2/3) /
  (2*3^(1/3)*Pi*r):

gg99_f := x -> gg99_f0(gg99_r(x)):

f := (rs, zeta, xt, xs0, xs1) ->
  gga_exchange(gg99_f, rs, zeta, xs0, xs1):
