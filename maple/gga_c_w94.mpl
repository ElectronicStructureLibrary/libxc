(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: work_gga_c *)

a := -1:
b := 11.8:
c :=  0.150670:
d := 11.02e-3/RS_FACTOR:

(* I added an abs() to this equation, as otherwise the f_num is
   complex for negative values of z. Of course it is not clear at all
   what was the original intential of Wilson, or if he even considered
   this problem. *)
f_num := z -> a*sqrt(1 - m_abs(1*z)^(5/3)):
f_den := (rs, xt) -> b + c*xt^(51/16) + d*xt^2*rs + rs:

(* Equation (25) *)
fa := (rs, z, xt, xs0, xs1) ->
  f_num(z)/f_den(rs, xt):

f := (rs, z, xt, xs0, xs1) ->
  f_num(z)/f_den(rs, xt):
