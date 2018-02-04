(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: work_gga_c *)

$include "gga_c_tca.mpl"

msinc := x -> piecewise(x = 0, 1, sin(x)/x):
aa    := Pi*(9*Pi/4)^(1/3):

fD := (rs, z, s) -> 1 - z^4*(1 - msinc(aa*s/rs)^2):

f := (rs, z, xt, xs0, xs1) -> 
  f_tcs(rs, z, xt)*fD(rs, z, X2S*2^(1/3)*xt):