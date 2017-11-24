(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: work_gga_x *)

$define gga_x_rpw86_params
$include "gga_x_pw86.mpl"

malpha := 0.02178:
mbeta  := 1.15:
muLV   := 0.8491/9:

f0 := s -> 
   + (1 + muLV*s^2)/(1 + malpha*s^6) 
   + malpha*s^6*f0_pw86(s)/(mbeta + malpha*s^6):

f  := x -> f0(X2S*x):

