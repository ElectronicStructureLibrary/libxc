(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: work_gga_c *)

$define lda_c_pw_params
$define lda_c_pw_modified_params
$include "lda_c_pw.mpl"

malpha := 2.804:
mgamma := 0.8098:

XX := s -> 1/(1 + malpha*s^2):
ff := s -> XX(s) + mgamma*(1 - XX(s)):

f := (rs, z, xt, xs0, xs1) -> f_pw(rs, z)*(
  + (1 + z)/2 * ff(X2S*xs0)
  + (1 - z)/2 * ff(X2S*xs1)
):