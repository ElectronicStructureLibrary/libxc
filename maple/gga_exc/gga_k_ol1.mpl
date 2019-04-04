(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)

ol1_c4 := 0.00677:

ol1_f0 := s -> 1 + (s^2/72 + ol1_c4*s)/K_FACTOR_C:
ol1_f  := x -> ol1_f0(x/2^(1/3)):

f := (rs, zeta, xt, xs0, xs1) -> gga_kinetic(ol1_f, rs, zeta, xs0, xs1):
