(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: work_gga_x *)

mkappa := KAPPA_PBE:
mmu    := 0.00361218645365094697:
malpha := 0.52:

f  := x -> 1 + mkappa*(1 - (1 + mmu*x^2/(malpha*mkappa))^(-malpha)):
