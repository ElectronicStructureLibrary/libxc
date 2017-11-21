(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: work_gga_x *)

c4 := 0.00677:

f0 := s -> 1 + (s^2/72 + c4*s)/K_FACTOR_C:
f  := x -> f0(x/2^(1/3)):