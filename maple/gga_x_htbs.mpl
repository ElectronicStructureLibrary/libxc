(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: work_gga_x *)

s1 := 0.6:
s2 := 2.6:

(* The equations to solve in order to obtain the coeficients cc are
  G(s1) = 0
  G(s2) = 1
 G'(s1) = 0
 G'(s2) = 0
G''(s1) = 0
G''(s2) = 0
*)

cc0 :=  s1^3*(s1^2 - 5*s1*s2 + 10*s2^2)/(s1 - s2)^5:
cc1 := -30*s1^2*s2^2/(s1 - s2)^5:
cc2 :=  30*s1*s2*(s1 + s2)/(s1 - s2)^5:
cc3 := -10*(s1^2 + 4*s1*s2 + s2^2)/(s1 - s2)^5:
cc4 :=  15*(s1 + s2)/(s1 - s2)^5:
cc5 := -6/(s1 - s2)^5:

$define gga_x_rpbe_params
$include "gga_x_rpbe.mpl"
$include "gga_x_wc.mpl"

g := s -> cc0 + cc1*s + cc2*s^2 + cc3*s^3 + cc4*s^4 + cc5*s^5:

f0 := s -> convert(piecewise(
   s < s1, f0_wc(s),
   s > s2, f0_rpbe(s),
   g(s)*f0_rpbe(s) + (1 - g(s))*f0_wc(s)
), 'Heaviside'):

f  := x -> f0(X2S*x):
