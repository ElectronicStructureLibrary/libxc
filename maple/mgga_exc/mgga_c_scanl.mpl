(*
 Copyright (C) 2019 D. Mejia-Rodriguez

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)

$define mgga_k_pcopt_params
$include "mgga_k_pc07.mpl"
$include "mgga_c_scan.mpl"

f := (rs, z, xt, xs0, xs1, us0, us1, ts0, ts1) ->
  scan_f(rs, z, xt, xs0, xs1, K_FACTOR_C*pc07_f(xs0,us0), 
    K_FACTOR_C*pc07_f(xs1,us1)):
