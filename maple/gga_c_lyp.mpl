(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: work_gga_c *)
(* prefix:
  gga_c_lyp_params *params;

  assert(p->params != NULL);
  params = (gga_c_lyp_params * )(p->params);
*)

Cf := 3/10 * (3*Pi^2)^(2/3):

omega := rs -> params_a_B*exp(-params_a_c*rs)/(1 + params_a_d*rs):
delta := rs -> (params_a_c + params_a_d/(1 + params_a_d*rs))*rs:

aux6 := 1/2^(8/3):
aux4 := aux6/4:
aux5 := aux4/(9*2):

t1 := (rs, z) ->
  -(1 - z^2)/(1 + params_a_d*rs):
t2 := (rs, z, xt) ->
  -xt^2*((1 - z^2)*(47 - 7*delta(rs))/(4*18) - 2/3):
t3 := (z) ->
  -Cf/2*(1 - z^2)*((1 + z)^(8/3) + (1 - z)^(8/3)):
t4 := (rs, z, xs0, xs1) ->
  aux4*(1 - z^2)*(5/2 - delta(rs)/18)*(xs0^2*(1 + z)^(8/3) + xs1^2*(1 - z)^(8/3)):
t5 := (rs, z, xs0, xs1) ->
  aux5*(1 - z^2)*(delta(rs) - 11)*(xs0^2*(1 + z)^(11/3) + xs1^2*(1 - z)^(11/3)):
t6 := (z, xs0, xs1) ->
  -aux6*(2/3*(xs0^2*(1 + z)^(8/3) + xs1^2*(1 - z)^(8/3))
  -(1 + z)^2*xs1^2*(1 - z)^(8/3)/4 - (1 - z)^2*xs0^2*(1 + z)^(8/3)/4):

f_lyp := (rs, z, xt, xs0, xs1) -> params_a_A*(t1(rs/RS_FACTOR, z) + omega(rs/RS_FACTOR)*(
  + t2(rs/RS_FACTOR, z, xt) + t3(z) + t4(rs/RS_FACTOR, z, xs0, xs1)
  + t5(rs/RS_FACTOR, z, xs0, xs1) + t6(z, xs0, xs1)
)):

f  := (rs, z, xt, xs0, xs1) -> f_lyp(rs, z, xt, xs0, xs1):

