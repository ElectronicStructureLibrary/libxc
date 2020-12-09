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

ml1_C := 6.187335:
ml1_b := [2.763169, 1.757515, 1.741397, 0.568985, 1.572202, 1.885389]:

ml1_alpha := z -> params_a_fc*((1 + z)^params_a_q + (1 - z)^params_a_q):
ml1_beta  := z -> (1 - z^2)^(1/3)/((1 + z)^(1/3) + (1 - z)^(1/3)):

z_thr2 := z -> my_piecewise5(
  1 + z <= p_a_zeta_threshold, p_a_zeta_threshold - 1,
  1 - z <= p_a_zeta_threshold, p_a_zeta_threshold + 1,
  z):

ml1_k := (rs, z) -> ml1_C*RS_FACTOR/rs * ml1_alpha(z)*ml1_beta(z):
   
ml1_Q := (rs, z) ->
  - ml1_b[1]/(1 + ml1_b[2]*ml1_k(rs, z))
  + ml1_b[3]*log(1 + ml1_b[4]/ml1_k(rs, z))/ml1_k(rs, z)
  + ml1_b[5]/ml1_k(rs, z)
  - ml1_b[6]/ml1_k(rs, z)^2:

(* screen for small spin densities to avoid divergences in the
  potentials.  Note that beta is zero for any polarized density and
  the whole expression for alpha*beta is symmetric in z.  Note also
  that in the expression for Q one divides by k that is zero for
  ferromagnetic densities. *)
  
ml1_f := (rs, z) -> 1/2*(RS_FACTOR/rs)^3 *
  my_piecewise3(1 - abs(z) <= p_a_zeta_threshold, 0, (1 - z^2)/4 * ml1_Q(rs, z_thr2(z))):

f := (rs, z) -> ml1_f(rs, z):
