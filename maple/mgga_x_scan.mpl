(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: work_mgga_x *)

p     := x -> X2S^2*x^2:
alpha := (x, t) -> (t - x^2/8)/K_FACTOR_C:

c1 := 0.667:
c2 := 0.8:
d  := 1.24:
f_alpha := a -> piecewise(a <= 1, exp(-c1*a/(1 - a)), -d*exp(c2/(1 - a))):

k1  := 0.065:
h1x := x -> 1 + k1*(1 - k1/(k1 + x)):

b2 := sqrt(5913/405000):
b1 := (511/13500)/(2*b2):
b3 := 1/2:
b4 := MU_GE^2/k1 - 1606/18225 - b1^2:
y  := (x, a) -> MU_GE*p(x) + b4*p(x)^2*exp(-b4*p(x)/MU_GE) + (b1*p(x) + b2*(1 - a)*exp(-b3*(1 - a)^2))^2:

a1 := 4.9479:
gx := x -> 1 - exp(-a1/sqrt(X2S*x)):

h0x := 1.174:
f   := (rs, x, t, u) -> (h1x(y(x, alpha(x, t)))*(1 - f_alpha(alpha(x, t))) + h0x*f_alpha(alpha(x, t)))*gx(x):
