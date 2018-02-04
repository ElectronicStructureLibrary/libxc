(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: work_gga_x *)

f0 := s -> 1 + B1*s*log(1 + s) + B2*s*log(1 + log(1 + s)):
f  := x -> f0(X2S*x):
