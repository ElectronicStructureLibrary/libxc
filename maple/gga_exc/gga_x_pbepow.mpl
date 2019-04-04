(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)

pbepow_kappa := KAPPA_PBE:
pbepow_mu    := 0.2195149727645171:
pbepow_m     := 100:

pbepow_gamm  := pbepow_m*pbepow_mu/pbepow_kappa:
pbepow_Cx    := pbepow_kappa/pbepow_m:

pbepow_f0 := s -> 1 + add(pbepow_Cx * (pbepow_gamm*s^2/(1 + pbepow_gamm*s^2))^i, i=1..pbepow_m):

pbepow_f  := x -> pbepow_f0(X2S*x):

f := (rs, z, xt, xs0, xs1) -> gga_exchange(pbepow_f, rs, z, xs0, xs1):
