(*
 Copyright (C) 2017 M.A.L. Marques
 Copyright (C) 2018 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: work_gga_x *)
(* prefix:
  gga_x_lspbe_params *params;
 
  assert(p->params != NULL);
  params = (gga_x_lspbe_params * )(p->params);
*)

f0_lspbe := s -> 1 + params_a_kappa*(1 - params_a_kappa/(params_a_kappa + params_a_mu*s^2)) - (params_a_kappa+1)*(1-exp(-params_a_alpha*s^2)):
f_lspbe  := x -> f0_lspbe(X2S*x):

f  := x -> f_lspbe(x):
