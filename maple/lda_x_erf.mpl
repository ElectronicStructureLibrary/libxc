(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: work_lda *)

$define lda_x_params
$include "lda_x.mpl"

(* attenuation function *)
aux1_erf := a -> sqrt(Pi)*erf(1/(2*a)):
aux2_erf := a -> exp(-1/(4*a^2)) - 1:
aux3_erf := a -> 2*a^2*aux2_erf(a) + 1/2:
attenuation_erf := a -> 1 - 8/3*a*(aux1_erf(a) + 2*a*(aux2_erf(a) - aux3_erf(a))):

a_cnst := (4/(9*Pi))^(1/3)*p_a_cam_omega/2:

lda_x_erf_s := (rs, z) -> convert(piecewise(z = -1, 0,
  (1 + z)^(4/3)*attenuation_erf(a_cnst*rs/(1 + z)^(1/3))
), 'Heaviside'):

f_lda_x_erf := (rs, z) -> lda_x_ax*(
   lda_x_erf_s(rs, z) + lda_x_erf_s(rs, -z)
)/rs:

f := (rs, z) -> f_lda_x_erf(rs, z):
