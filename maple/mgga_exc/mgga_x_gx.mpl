(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: work_mgga_x *)

malpha := (x, t) -> (t - x^2/8)/K_FACTOR_C:

gx_cx0 := 4/3*(2/Pi)^(1/3):
gx_cx1 := X_FACTOR_C:

gx_c0 :=  0.827411:
gx_c1 := -0.643560:

gx_gx0 := a ->
  + gx_cx0/gx_cx1
  + a*(gx_c0 + gx_c1*a)/(1.0 + (gx_c0 + gx_c1 - 1)*a) * (1 - gx_cx0/gx_cx1):

gx_alphainf := 0.852:
gx_gx1 := a ->
  1 + (1 - gx_alphainf)*(1 - a)/(1 + a):

f_gx_a := a->
  + gx_gx0(a)*Heaviside(1 - a)
  + gx_gx1(a)*Heaviside(a - 1):

f := (rs, x, t, u) ->
  f_gx_a(malpha(x, t)):
