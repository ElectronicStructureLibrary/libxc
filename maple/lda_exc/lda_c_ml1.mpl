(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: lda_exc *)
(* prefix:
  lda_c_ml1_params *params;

  assert(p->params != NULL);
  params = (lda_c_ml1_params * )(p->params);
*)

ml1_CC := 6.187335:
ml1_bb := [2.763169, 1.757515, 1.741397, 0.568985, 1.572202, 1.885389]:

ml1_alpha := z -> params_a_fc*((1 + z)^params_a_q + (1 - z)^params_a_q):
ml1_beta  := z -> (1 - z^2)^(1/3)/((1 + z)^(1/3) + (1 - z)^(1/3)):

(* screen for small spin densities to avoid divergences in the
  potentials.  Note that beta is zero for any polarized density and
  the expression for alpha*beta is symmetric in z *)
  
z_thr2 := z -> my_piecewise5(
  1 + z <= p_a_zeta_threshold, p_a_zeta_threshold - 1,
  1 - z <= p_a_zeta_threshold, p_a_zeta_threshold + 1,
  z):

ml1_kk := (rs, z) -> ml1_CC*RS_FACTOR/rs *
   my_piecewise3(1 - m_abs(z) <= p_a_zeta_threshold, 0, ml1_alpha(z_thr2(z))*ml1_beta(z_thr2(z))):
   
ml1_QQ := (rs, z) ->
  - ml1_bb[1]/(1 + ml1_bb[2]*ml1_kk(rs, z))
  + ml1_bb[3]*log(1 + ml1_bb[4]/ml1_kk(rs, z))/ml1_kk(rs, z)
  + ml1_bb[5]/ml1_kk(rs, z)
  - ml1_bb[6]/ml1_kk(rs, z)^2:

ml1_f := (rs, z) -> 1/2*(RS_FACTOR/rs)^3 * (1 - z^2)/4 * ml1_QQ(rs, z):

f := (rs, z) -> ml1_f(rs, z):
