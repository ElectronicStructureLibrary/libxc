(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: lda_exc *)
(* prefix:
  lda_xc_ksdt_params *params;

  assert(p->params != NULL);
  params = (lda_xc_ksdt_params * )(p->params);
*)

ksdt_g := [2/3, -0.0139261, 0.183208]:
ksdt_l := [1.064009, 0.572565]:

ksdt_alpha := (t, rs) ->
      2 - (ksdt_g[1] + ksdt_g[2]*rs)/(1 + ksdt_g[3]*rs)*exp(-t*(ksdt_l[1] + ksdt_l[2]*t*sqrt(rs))):

ksdt_phi := (z, malpha) ->
      (opz_pow_n(z, malpha) + opz_pow_n(-z, malpha) - 2)/(2^malpha - 2):


ksdt_lambda := (4/(9*Pi))^(1/3):
ksdt_pa0    := 1/(Pi*ksdt_lambda):
ksdt_pa     := [0.750, 3.043630, -0.0922700, 1.703500, 8.310510, 5.11050]:

ksdt_faa := t -> ksdt_pa0*tanh(1/t)*(ksdt_pa[1] + ksdt_pa[2]*t^2 + ksdt_pa[3]*t^3 + ksdt_pa[4]*t^4)
    / (1 + ksdt_pa[5]*t^2 + ksdt_pa[6]*t^4):
ksdt_fbb := (b, t) -> tanh(1/sqrt(t))*(b[1] + b[2]*t^2 + b[3]*t^4)/(1 + b[4]*t^2 + b[5]*t^4):
ksdt_fdd := (d, t) -> ksdt_fbb(d, t):
ksdt_fee := (e, t) -> tanh(1/t)*(e[1] + e[2]*t^2 + e[3]*t^4)/(1 + e[4]*t^2 + e[5]*t^4):
ksdt_fcc := (c, e, t) -> (c[1] + c[2]*exp(-c[3]/t))*ksdt_fee(e, t):

ksdt_fxc := (omega, b, c, d, e, rs, t) ->
    -(omega*ksdt_faa(t) + ksdt_fbb(b, t)*sqrt(rs) + ksdt_fcc(c, e, t)*rs)/(rs*(1 + ksdt_fdd(d, t)*sqrt(rs) + ksdt_fee(e, t)*rs)):

(* (T/T_F)*opz_pow_n(z, 2/3) *)
ksdt_mtt := (rs, z) ->
    2*(4/(9*Pi))^(2/3)*params_a_ksdt_T*rs^2*(1 + params_a_ksdt_thetaParam*z)^(2/3):

f := (rs, z) ->
  + ksdt_fxc(1,
        params_a_ksdt_pb_0_, params_a_ksdt_pc_0_, params_a_ksdt_pd_0_, params_a_ksdt_pe_0_,
        rs, ksdt_mtt(rs, z))*(1 - ksdt_phi(z, ksdt_alpha(ksdt_mtt(rs, z), rs)))
  + ksdt_fxc(2^(1/3),
        params_a_ksdt_pb_1_, params_a_ksdt_pc_1_, params_a_ksdt_pd_1_, params_a_ksdt_pe_1_,
        rs, ksdt_mtt(rs, z)/2^(2/3))*ksdt_phi(z, ksdt_alpha(ksdt_mtt(rs, z), rs)):
