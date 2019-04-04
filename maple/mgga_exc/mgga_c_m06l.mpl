(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: work_mgga_c *)
(* prefix:
  mgga_c_m06l_params *params;

  assert(p->params != NULL);
  params = (mgga_c_m06l_params * )(p->params);
*)

$include "mgga_c_vsxc.mpl"
$include "mgga_c_m05.mpl"

f_m06 := (rs, z, xs0, xs1, ts0, ts1) ->
  + f_m05(rs, z, xs0, xs1, ts0, ts1)
  + f_vsxc(rs, z, xs0, xs1, ts0, ts1):

f := (rs, z, xt, xs0, xs1, ts0, ts1, us0, us1) ->
  f_m06(rs, z, xs0, xs1, ts0, ts1):

