(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)
(* prefix:
  gga_x_sfat_pbe_params *params;

  assert(p->params != NULL);
  params = (gga_x_sfat_pbe_params * )(p->params);
*)

$define gga_x_pbe_params
$include "gga_x_pbe.mpl"
$include "gga_x_sfat.mpl"

ityh_enhancement := xs -> pbe_f(xs):
