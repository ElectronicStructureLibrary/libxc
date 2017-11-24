(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: work_lda *)

C1 := -0.0603:
C2 :=  0.0175:
C3 := -0.00053:

f := (rs, zeta) -> C1 + C2*n_total(rs)^(-1/3) + C3*n_total(rs)^(-2/3):

