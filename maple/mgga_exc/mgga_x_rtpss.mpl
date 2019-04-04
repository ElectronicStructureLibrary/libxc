(*
 Copyright (C) 2017 M.A.L. Marques
 Copyright (C) 2018 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: work_mgga_x *)
(* prefix:
  mgga_x_rtpss_params *params;
 
  assert(pt->params != NULL);
  params = (mgga_x_rtpss_params * ) (pt->params);
*)

(* These are used within the tpss_x routine *)
ff     := z -> 2:
mkappa := (x, t) -> params_a_kappa:

$include "tpss_x.mpl"

(* Equation (6) *)

f   := (rs, x, t, u) -> 1 + mkappa(x, t)*(1 - exp(-fx(x, t)/mkappa(x,t))):
