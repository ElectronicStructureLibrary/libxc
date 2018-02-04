(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: work_gga_x *)

c     :=  0.7168:
alpha :=  2.804:
d     := 28.23705740248932030511071641312341561894: (* POW(CBRT(4/3) * 2*M_PI/3, 4) *)

csi  := s -> (3/2 * LambertW(s^(3/2) / (2*sqrt(6))))^(2/3):
fb   := s -> Pi/3 * s/(csi(s) * (d + csi(s)^2)^(1/4)):
flaa := s -> (1 + c*s^2)/(1 + c*s^2/fb(s)):
XX   := s -> 1 - alpha*s^2/(1 + alpha*s^2):

f    := x ->  XX(X2S*x) + (1 - XX(X2S*x))*flaa(X2S*x):
