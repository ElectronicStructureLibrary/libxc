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


f := (rs, z, xt, xs0, xs1, u0, u1, t0, t1) -> 
  mgga_exchange(scan_f, rs, z, xs0, xs1, u0, u1, 
    K_FACTOR_C*pc07_f(xs0,u0), K_FACTOR_C*pc07_f(xs1,u1)):
