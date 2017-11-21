(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: work_mgga_x *)

a2 := 146/2025:
a3 := -73/405:
a4 := 0.131957187845257783631757384393: (* DD + 100/(81*81*kappa) *)

kappa := 0.804:

xx := (p, qt) -> MU_GE*p + a2*qt*qt + a3*qt*p + a4*p*p:

f := (rs, x, t, u) -> 1 + kappa - \
  kappa^2/(kappa + xx(X2S^2*x^2, 6*X2S^2*t - 9/20 - X2S^2*x^2/12)):
