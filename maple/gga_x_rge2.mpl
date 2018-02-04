(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: work_gga_x *)

kappa := 0.8040:

den := s -> kappa + 1*MU_GE*s^2 + MU_GE^2*s^4/kappa:
f0  := s -> 1 + kappa * (1 - kappa/den(s)):
f   := x -> f0(X2S * x):
