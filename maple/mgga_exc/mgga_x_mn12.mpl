(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: work_mgga_c *)
(* prefix:
  mgga_x_mn12_params *params;

  assert(p->params != NULL);
  params = (mgga_x_mn12_params * ) (p->params);
*)

$define lda_x_params
$include "lda_x.mpl"

omega_x := 2.5:
gamma_x := 0.004:

vx := (rs, z) -> 1/(1 + rs/(omega_x*RS_FACTOR)*(2/(1 + z))^(1/3)):
ux := x -> gamma_x*x^2/(1 + gamma_x*x^2):
wx := t -> (K_FACTOR_C - t)/(K_FACTOR_C + t):

pol1  := t-> params_a_c[ 1] + params_a_c[ 2]*wx(t) + params_a_c[ 3]*wx(t)^2 + params_a_c[ 4]*wx(t)^3
  + params_a_c[ 5]*wx(t)^4 + params_a_c[ 6]*wx(t)^5:
pol2  := t-> params_a_c[ 7] + params_a_c[ 8]*wx(t) + params_a_c[ 9]*wx(t)^2 + params_a_c[10]*wx(t)^3
   + params_a_c[11]*wx(t)^4:
pol3  := t-> params_a_c[12] + params_a_c[13]*wx(t) + params_a_c[14]*wx(t)^2 + params_a_c[15]*wx(t)^3:
pol4  := t-> params_a_c[16] + params_a_c[17]*wx(t) + params_a_c[18]*wx(t)^2:
pol5  := t-> params_a_c[19] + params_a_c[20]*wx(t) + params_a_c[21]*wx(t)^2 + params_a_c[22]*wx(t)^3
   + params_a_c[23]*wx(t)^4:
pol6  := t-> params_a_c[24] + params_a_c[25]*wx(t) + params_a_c[26]*wx(t)^2 + params_a_c[27]*wx(t)^3:
pol7  := t-> params_a_c[28] + params_a_c[29]*wx(t) + params_a_c[30]*wx(t)^2:
pol8  := t-> params_a_c[31] + params_a_c[32]*wx(t) + params_a_c[33]*wx(t)^2 + params_a_c[34]*wx(t)^3:
pol9  := t-> params_a_c[35] + params_a_c[36]*wx(t) + params_a_c[37]*wx(t)^2:
pol10 := t-> params_a_c[38] + params_a_c[39]*wx(t) + params_a_c[40]*wx(t)^2:

FMN12 := (rs, z, x, t) ->
  + pol1(t)
  + pol2(t)*ux(x)
  + pol3(t)*ux(x)^2
  + pol4(t)*ux(x)^3
  + pol5(t)*vx(rs, z)
  + pol6(t)*ux(x)*vx(rs, z)
  + pol7(t)*ux(x)^2*vx(rs, z)
  + pol8(t)*vx(rs, z)^2
  + pol9(t)*ux(x)*vx(rs, z)^2
  + pol10(t)*vx(rs, z)^3:

f_spin := (rs, z, x, t) ->
  lda_x_ax*(1 + z)^(4/3)/rs*FMN12(rs, z, x, t):

f_mn12 := (rs, z, xt, xs0, xs1, ts0, ts1) ->
  f_spin(rs, z, xs0, ts0) + f_spin(rs, -z, xs1, ts1):

f := (rs, z, xt, xs0, xs1, ts0, ts1, us0, us1) ->
  f_mn12(rs, z, xt, xs0, xs1, ts0, ts1):
