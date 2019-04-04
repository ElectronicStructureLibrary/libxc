(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)

cap_Ax       := -3/4*(3/Pi)^(1/3):
cap_mu       := 0.2195149727645171:
cap_alpha    := -cap_Ax*cap_mu:
cap_c        := cap_alpha/(3*Pi^2)^(1/3):
cap_alphaoAx := -cap_mu:

cap_f0 := s -> 1 - cap_alphaoAx*s*log(1 + s)/(1 + cap_c*log(1 + s)):
cap_f  := x -> cap_f0(X2S*x):

f := (rs, z, xt, xs0, xs1) -> gga_exchange(cap_f, rs, z, xs0, xs1):

