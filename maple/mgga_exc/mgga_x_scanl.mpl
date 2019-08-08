(*
 Copyright (C) 2019 D. Mejia-Rodriguez

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)
(* prefix:
  mgga_x_scan_params *params;

  assert(p->params != NULL);
  params = (mgga_x_scan_params * )(p->params);
*)

$define mgga_k_pcopt_params
$include "mgga_k_pc07.mpl"
$include "mgga_x_scan.mpl"

scanl_alpha := (x, u) -> pc07_f(x, u) - pc07_f_W(x):

scanl_f   := (x, u) -> (scan_h1x(scan_y(x, scanl_alpha(x, u)))*(1 - scan_f_alpha(scanl_alpha(x, u)))
  + scan_h0x*scan_f_alpha(scanl_alpha(x, u)))*scan_gx(x):

f := (rs, z, xt, xs0, xs1, u0, u1, t0, t1) -> mgga_exchange(scanl_f, rs, z, xs0, xs1, u0, u1, t0, t1):
