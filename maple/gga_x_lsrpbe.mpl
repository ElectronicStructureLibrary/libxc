(*
 Copyright (C) 2017 M.A.L. Marques
 Copyright (C) 2018 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: work_gga_x *)
(* prefix:
  gga_x_lsrpbe_params *params;

  assert(p->params != NULL);
  params = (gga_x_lsrpbe_params * )(p->params);
*)

f0_lsrpbe := s -> 1 + params_a_kappa * (
  1 - exp(-params_a_mu*s^2/params_a_kappa)
) - (params_a_kappa+1)*(1 - exp(-params_a_alpha*s^2)):
f_lsrpbe  := x -> f0_lsrpbe(X2S*x):

f       := x -> f_lsrpbe(x):
