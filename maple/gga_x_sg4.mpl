(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: work_gga_x *)

mu1     := 0.042:
mu2     := 0.26 - mu1:
nu_MGE4 := -0.195:
k2      := -mu2^2/nu_MGE4:
k1      := 0.804 - k2:

f0 := s -> 1 + k1 + k2 
   - k1*(1 - mu1*s^2/k1)/(1 - (mu1*s^2/k1)^5)
   - k2/(1 + mu2*s^2/k2):

f:= x -> f0(X2S*x): 