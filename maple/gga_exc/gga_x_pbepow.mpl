(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)

pbepow_kappa := KAPPA_PBE:
pbepow_mu    := 0.2195149727645171:

(* pbepow_m has to be written in float format otherwise the code is converted to C
   in a numerically unstable way *)
pbepow_m     := 100.0:

pbepow_gamm  := pbepow_mu/pbepow_kappa:

(* The geometrical series 1 + pbepow_kappa/pbepow_m * add(pbepow_yy(X2S*s)^i, i=1..pbepow_m)
   is summed analytically and simplified *)
pbepow_f0 := s -> 1 - pbepow_kappa/pbepow_m
  + pbepow_kappa*(1/pbepow_m + pbepow_gamm*s^2)*(1 - (pbepow_gamm*s^2/(1/pbepow_m + pbepow_gamm*s^2))^(pbepow_m + 1)):

pbepow_f  := x -> pbepow_f0(X2S*x):

f := (rs, z, xt, xs0, xs1) -> gga_exchange(pbepow_f, rs, z, xs0, xs1):
