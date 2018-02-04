(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: work_gga_x *)

lambda := y -> 1/2*(1 + (1 - y^2)*log((1 + y)/abs(1 - y))/(2*y)):

f := x -> 1 + lambda(X2S*x/6)*x^2/(8*K_FACTOR_C):