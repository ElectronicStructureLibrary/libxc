(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: work_gga_x *)

mkappa := 0.4604:
mmu    := 0.354546875:

f0 := s -> 1 + mkappa*(1 - mkappa/(mkappa + mmu*s^2)):
f  := x -> f0(X2S_2D*x):
