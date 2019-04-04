(*
 Copyright (C) 2017 M.A.L. Marques
 Copyright (C) 2018 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: work_mgga_x *)
(* prefix:
  mgga_x_msb_params *params;

  assert(pt->params != NULL);
  params = (mgga_x_msb_params * ) (pt->params);
*)

$include "mgga_x_ms.mpl"

beta := (t,x) -> alpha(t,x)*K_FACTOR_C/(t+K_FACTOR_C):

f := (rs, x, t, u) -> f0(X2S^2*x^2, 0) + \
  fa(beta(t,x))*(f0(X2S^2*x^2, params_a_c) - f0(X2S^2*x^2, 0)):
