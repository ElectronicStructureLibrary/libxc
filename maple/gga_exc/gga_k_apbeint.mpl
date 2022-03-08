(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)
(* prefix:
  assert(p->params != NULL);
  const gga_k_apbeint_params * const params = (gga_k_apbeint_params * const)(p->params);
*)

$include "gga_x_pbeint.mpl"

f := (rs, z, xt, xs0, xs1) -> gga_kinetic(pbeint_f, rs, z, xs0, xs1):
