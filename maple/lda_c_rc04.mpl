(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: work_lda *)

AA := -0.655868:
BB :=  4.888270:
CC :=  3.177037:
DD :=  0.897889:

phi := z -> 1/2*((1 + z)^(2/3) + (1 - z)^(2/3)):

f_rc04 := (rs, zeta) -> phi(zeta)^3 * (AA*arctan(BB + CC*rs) + DD)/rs:

f      := (rs, zeta) -> f_rc04(rs, zeta):