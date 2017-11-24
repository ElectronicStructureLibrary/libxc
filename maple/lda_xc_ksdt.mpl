(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: work_lda *)
(* prefix:
  lda_xc_ksdt_params *params;

  assert(p->params != NULL);
  params = (lda_xc_ksdt_params * )(p->params);
*)

g := [2/3, -0.0139261, 0.183208]:
l := [1.064009, 0.572565]:

alpha := (t, rs) ->
      2 - (g[1] + g[2]*rs)/(1 + g[3]*rs)*exp(-t*(l[1] + l[2]*t*sqrt(rs))):


phi := (malpha, z) -> ((1 + z)^malpha + (1 - z)^malpha - 2)/(2^malpha - 2):


lambda := (4/(9*Pi))^(1/3):
a0     := 1/(Pi*lambda):

a     := [0.750, 3.043630, -0.0922700, 1.703500, 8.310510, 5.11050]:

aa := t -> a0*tanh(1/t)*(a[1] + a[2]*t^2 + a[3]*t^3 + a[4]*t^4)/(1 + a[5]*t^2 + a[6]*t^4):
bb := (b, t) -> tanh(1/sqrt(t))*(b[1] + b[2]*t^2 + b[3]*t^4)/(1 + b[4]*t^2 + b[5]*t^4):
dd := (d, t) -> bb(d, t):
ee := (e, t) -> tanh(1/t)*(e[1] + e[2]*t^2 + e[3]*t^4)/(1 + e[4]*t^2 + e[5]*t^4):
cc := (c, e, t) -> (c[1] + c[2]*exp(-c[3]/t))*ee(e, t):

fxc := (omega, b, c, d, e, rs, t) ->
    -(omega*aa(t) + bb(b, t)*sqrt(rs) + cc(c, e, t)*rs)/(rs*(1 + dd(d, t)*sqrt(rs) + ee(e, t)*rs)):

temp := max(params_a_T, 1e-8):
tt := rs -> 2*(4/(9*Pi))^(2/3)*temp*rs^2:

f := (rs, z) ->
  + fxc(1,
        params_a_b_0_, params_a_c_0_, params_a_d_0_, params_a_e_0_,
        rs, tt(rs))*(1 - phi(alpha(tt(rs), rs), z))
  + fxc(2^(1/3),
        params_a_b_1_, params_a_c_1_, params_a_d_1_, params_a_e_1_,
        rs, tt(rs))*phi(alpha(tt(rs), rs), z):
