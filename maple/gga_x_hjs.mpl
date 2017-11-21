(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: work_gga_c *)
(* prefix:
  gga_x_hjs_params *params;

  assert(p->params != NULL);
  params = (gga_x_hjs_params * )(p->params);
*)

AA :=  0.757211:
BB := -0.106364:
CC := -0.118649:
DD :=  0.609650:

fH := s -> add(params_a_a[i]*s^(1+i), i=1..6)/(1 + add(params_a_b[i]*s^i, i=1..9)):

mzeta   := s -> s^2*fH(s):
meta    := s -> AA + mzeta(s):
mlambda := s -> DD + mzeta(s):
mchi    := (rs, z, s) -> nu(rs, z)/sqrt(mlambda(s) + nu(rs, z)^2):

fF := (rs, z, s) ->
  1 - s^2/(27*CC*(1 + s^2/4)) - mzeta(s)/(2*CC):

fG := (rs, z, s) ->
  - 2/5  * CC*fF(rs, z, s)*mlambda(s)
  - 4/15 * BB*mlambda(s)^2
  - 6/5  * AA*mlambda(s)^3
  - mlambda(s)^(7/2)*(4/5*sqrt(Pi) + 12/5*(sqrt(mzeta(s)) - sqrt(meta(s)))):

f1 := (rs, z, s) ->
   + AA
   - 4/9 * BB*(1 - mchi(rs, z, s))/mlambda(s)
   - 2/9 * CC*fF(rs, z, s)*(2 - 3*mchi(rs, z, s) + mchi(rs, z, s)^3)/mlambda(s)^2
   - 1/9 * fG(rs, z, s)*(8 - 15*mchi(rs, z, s) + 10*mchi(rs, z, s)^3 - 3*mchi(rs, z, s)^5)/mlambda(s)^3
   + 2*nu(rs, z)*(sqrt(mzeta(s) + nu(rs, z)^2) -  sqrt(meta(s) + nu(rs, z)^2))
   + 2*mzeta(s)*log((nu(rs, z) + sqrt(mzeta(s) + nu(rs, z)^2))/(nu(rs, z) + sqrt(mlambda(s) + nu(rs, z)^2)))
   - 2*meta(s)*log((nu(rs, z) + sqrt(meta(s) + nu(rs, z)^2))/(nu(rs, z) + sqrt(mlambda(s) + nu(rs, z)^2))):

f_lda := (rs, z) -> -X_FACTOR_C*RS_FACTOR*((1 + z)/2)^(4/3)/rs:

f := (rs, z, xt, xs0, xs1) ->
  f_lda(rs, z)*f1(rs, z, X2S*xs0) + f_lda(rs, -z)*f1(rs, -z, X2S*xs1):
