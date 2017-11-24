(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: work_gga_c *)
(* prefix:
  gga_x_n12_params *params;

  assert(p->params != NULL);
  params = (gga_x_n12_params * )(p->params);
*)

omega_x := 2.5:
gamma_x := 0.004:

rss := (rs, z) -> rs * (2/(1 + z))^(1/3):

vx := rs -> 1/(1 + (1/(RS_FACTOR*omega_x))*rs):
ux := x -> gamma_x*x^2/(1 + gamma_x*x^2):

FN12 := (rs, x) -> 
  + add(1*params_a_CC_0_[i+1]*ux(x)^i, i=0..3)
  + add(1*params_a_CC_1_[i+1]*ux(x)^i, i=0..3) * vx(rs)
  + add(1*params_a_CC_2_[i+1]*ux(x)^i, i=0..3) * vx(rs)^2
  + add(1*params_a_CC_3_[i+1]*ux(x)^i, i=0..3) * vx(rs)^3:

f_n12 := (rs, z, xt, xs0, xs1) -> -X_FACTOR_C*RS_FACTOR*(
  + (1 + z)*FN12(rss(rs,  z), xs0)/(2*rss(rs,  z))
  + (1 - z)*FN12(rss(rs, -z), xs1)/(2*rss(rs, -z))
):

f  := (rs, z, xt, xs0, xs1) -> f_n12(rs, z, xt, xs0, xs1):