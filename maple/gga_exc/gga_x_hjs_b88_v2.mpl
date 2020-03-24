(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)

$include "gga_x_hjs.mpl"

params_a_a := [0.0253933, -0.0673075, 0.0891476, -0.0454168, -0.00765813, 0.0142506]:
params_a_b := [-2.6506, 3.91108, -3.31509, 1.54485, -0.198386,
  -0.136112, 0.0647862, 0.0159586, -0.000245066]:

hjs2_xi := 1/(exp(20) - 1):
hjs2_fs := s -> -log((exp(-s) + hjs2_xi)/(1 + hjs2_xi)):

hjs_fx := (rs, z, x) -> hjs_f1(rs, z, hjs2_fs(X2S*x)):
