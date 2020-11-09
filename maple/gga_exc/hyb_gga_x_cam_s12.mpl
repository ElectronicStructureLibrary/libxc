(*
 Copyright (C) 2020 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)
(* prefix:
  hyb_gga_x_cam_s12_params *params;

  assert(p->params != NULL);
  params = (hyb_gga_x_cam_s12_params * )(p->params);
*)

$include "gga_x_s12.mpl"
$include "gga_x_ityh.mpl"

ityh_enhancement := xs  -> s12g_f(xs):
