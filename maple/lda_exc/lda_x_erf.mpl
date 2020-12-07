(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: lda_exc *)

$include "attenuation.mpl"

a_cnst := (4/(9*Pi))^(1/3)*p_a_hyb_omega_0_/2:

xerf_ax := -RS_FACTOR*X_FACTOR_C/2^(4/3):
xerf_spin := (rs, z) ->
  xerf_ax*opz_pow_n(z, 4/3)/rs * attenuation_erf(a_cnst*rs/opz_pow_n(z, 1/3)):

xerf_f := (rs, z) -> xerf_spin(rs, z) + xerf_spin(rs, -z):

f := (rs, z) -> xerf_f(rs, z):
