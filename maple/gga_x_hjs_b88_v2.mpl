(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: work_gga_c *)
(* prefix:
  gga_x_hjs_b88_v2_params *params;
 
  assert(p->params != NULL);
  params = (gga_x_hjs_b88_v2_params * )(p->params);
*)


params_a_a := [0.0253933, -0.0673075, 0.0891476, -0.0454168, -0.00765813, 0.0142506]:
params_a_b := [-2.6506, 3.91108, -3.31509, 1.54485, -0.198386, 
  -0.136112, 0.0647862, 0.0159586, -0.000245066]:


xi := 1/(exp(20) - 1):
fs := s -> -log((exp(-s) + xi)/(1 + xi)):

$include "gga_x_hjs.mpl"

f := (rs, z, xt, xs0, xs1) -> 
  f_lda(rs, z)*f1(rs, z, fs(X2S*xs0)) + f_lda(rs, -z)*f1(rs, -z, fs(X2S*xs1)):