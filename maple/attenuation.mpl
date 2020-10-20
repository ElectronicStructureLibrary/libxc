(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* error function:
    Toulouse et al, Int. J. of Quant. Chem. 100, 1047 (2004); doi:10.1002/qua.20259
    Tawada et al, J. Chem. Phys. 120, 8425 (2004); doi:10.1063/1.1688752
*)
att_erf_aux1 := a -> sqrt(Pi)*erf(1/(2*a)):
att_erf_aux2 := a -> exp(-1/(4*a^2)) - 1:
att_erf_aux3 := a -> 2*a^2*att_erf_aux2(a) + 1/2:
attenuation_erf := a ->
  1 - 8/3*a*(att_erf_aux1(a) + 2*a*(att_erf_aux2(a) - att_erf_aux3(a))):

attenuation_erf_f2 := a ->
  1 + 24*a^2*((20*a^2 - 64*a^4)*exp(-1/(4*a^2)) - 3 - 36*a^2 + 64*a^4 + 10*sqrt(Pi)*a*erf(1/(2*a))):

attenuation_erf_f3 := a ->
  1 + 8*a/7*(
      + (-8*a + 256*a^3 - 576*a^5 + 3849*a^7 - 122880*a^9)*exp(-1/(4*a^2))
      + 24*a^3*(-35 + 224*a^2 - 1440*a^4 + 5120*a^6)
      + 2*sqrt(Pi)*(-2 + 60*a^2)*erf(1/(2*a))):

(* erf_gau - screening function = + 2 mu/sqrt(pi) exp(-mu^2 r^2)
    Song et al, J. Chem. Phys. 127, 154109 (2007); doi:10.1063/1.2790017
    You can recover the result in Int. J. of Quant. Chem. 100, 1047 (2004)
    by putting a = a/sqrt(3) and multiplying the whole attenuation function by -sqrt(3)
*)
att_gau_aux1 := (a) -> sqrt(Pi)*erf(1/(2*a)):
att_gau_aux2 := (a) -> exp(-1/(4*a^2)) - 1:
attenuation_gau := (a) ->
  -8/3*a*(att_gau_aux1(a) + 2*a*att_gau_aux2(a)*(1 - 8*a^2) - 4*a):

(* yukawa
    Akinaga and Ten-no, Chem. Phys. Lett. 462, 348 (2008); doi:10.1016/j.cplett.2008.07.103
*)
att_yuk_aux1 := a -> arctan(1, a):
att_yuk_aux2 := a -> log(1.0 + 1.0/a^2):
att_yuk_aux3 := a -> a^2 + 1:
attenuation_yukawa := a -> my_piecewise3(
    a > 50, 1/(9*a^2) - 1/(30*a^4),
    1.0 - 8/3*a*(att_yuk_aux1(a) + a/4*(1 - (att_yuk_aux3(a) + 2)*att_yuk_aux2(a)))
):

