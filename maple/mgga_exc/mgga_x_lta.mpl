(*
 Copyright (C) 2017 M.A.L. Marques
               2019 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)

lta_a1 := 2^(-1/5)*(10/(3*(6*Pi*Pi)^(2/3)))^(4/5): (* (2/C_F)^(4/5) *):

lta_f := (x, u, t) -> lta_a1*t^(4/5):

f := (rs, z, xt, xs0, xs1, u0, u1, t0, t1) -> mgga_exchange(lta_f, rs, z, xs0, xs1, u0, u1, t0, t1):
