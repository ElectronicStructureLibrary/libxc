(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: work_mgga_x *)
(* prefix:
  mgga_x_m06l_params *params;

  assert(pt->params != NULL);
  params = (mgga_x_m06l_params * )(pt->params);
*)

$define gga_x_pbe_params
$include "gga_x_pbe.mpl"
$include "gvt4.mpl"

malpha  := 0.00186726:
coeff_d := params_a_d:

(* there is a factor if 2 in the definition of z, as in Theor. Chem. Account 120, 215 (2008) *)
(* A MINUS was missing in Eq. (7) of the paper *)

f := (rs, x, t, u) ->
  + f_pbe(x)*mgga_series_w(params_a_a, 12, t)
  + gtv4(malpha, coeff_d, x, 2*(t - K_FACTOR_C)):
