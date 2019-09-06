(*
 Copyright (C) 2017 M.A.L. Marques
               2019 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)
$include "mgga_xc_b98.mpl"

(* Fitted parameters, given below Eq. 50 *)
edmgga_mu := 0.21:
edmgga_a  := 0.704:

(* Eq. 40: F_x'(0) *)
edmgga_fp0   := 27*edmgga_mu/50:
(* Eq. 27: lim_(x -> -infinity) F_x(x) *)
edmgga_fminf := 1/3 * (4*Pi^2/3)^(1/3):
(* Eq. 44 *)
edmgga_b := 2*Pi/9 * sqrt(6/5):

(* Eq. 47 *)
edmgga_c1 := edmgga_fminf:
(* Eq. 48 *)
edmgga_c2 := 1 - edmgga_c1:
(* Eq. 49 *)
edmgga_c3 := edmgga_c2*sqrt(2*edmgga_a)/(2*edmgga_b):

(* Eq. 50. There's some problem in the equation, since the first line,
1/edmgga_c3^3 * (edmgga_c3^2 - edmgga_fp0/(2*edmgga_b^2)) * edmgga_c2,
does not agree with the second line. The first line evaluates to a
negative value, which makes the functional blow up. The second line
yields values that are close to the reported ones, but Tao might have
used a more accurate numerical value in his implementation since the
agreement is not perfect. Unfortunately there's no way to check this
as he has passed away. *)
edmgga_c4 := 1/edmgga_c3^3 * (edmgga_c3^2 - 0.09834*edmgga_mu):

(* Eq. 46 *)
edmgga_x := Qb -> edmgga_a*Qb + sqrt(1 + (edmgga_a*Qb)^2):

(* Eq. 45 *)
edmgga_f_x :=  x -> edmgga_c1 + edmgga_c2*x/(1 + edmgga_c3*sqrt(x)*arcsinh(edmgga_c4*(x - 1))):

edmgga_f := (x, u, t) -> edmgga_f_x(edmgga_x(b98_q(x, u, t))):

f := (rs, z, xt, xs0, xs1, u0, u1, t0, t1) ->
  mgga_exchange(edmgga_f, rs, z, xs0, xs1, u0, u1, t0, t1):
