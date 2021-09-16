(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)
(* prefix:
  gga_c_apbe_bfit_params *params;

  assert(p->params != NULL);
  params = (gga_c_apbe_bfit_params * )(p->params);
*)

$include "gga_c_pbe.mpl"

params_a_gamma := (1 - log(2))/Pi^2: (* 0.031090690869654895034 *)
params_a_BB    := 1:

bfit_b0   := 3*params_a_bfit_mux/Pi^2:
bfit_binf := 3*(7/81)/Pi^2:
p_a_hyb_omega_0_ := params_a_bfit_omega:

bfit_beta := (mnu) ->
  (bfit_b0 + params_a_bfit_a1*mnu + params_a_bfit_a2*mnu^2 + bfit_binf*params_a_bfit_a3*mnu^3)/
  (1 + params_a_bfit_a1*mnu + params_a_bfit_a2*mnu^2 + params_a_bfit_a3*mnu^3):

(* this is the fitted version of beta, Eq. (19) *)
mbeta := (rs, t) -> bfit_beta(nu(rs, 1)):