(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: work_mgga_x *)
(* prefix:
  mgga_x_m05_params *params;

  assert(pt->params != NULL);
  params = (mgga_x_m05_params * )(pt->params);
*)

$define gga_x_pbe_params
$include "gga_x_pbe.mpl"

f := (rs, x, t, u) ->
  + params_a_csi_HF*f_pbe(x)*mgga_series_w(params_a_a, 12, t):
