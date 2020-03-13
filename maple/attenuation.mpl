(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* error function:
    Int. J. of Quant. Chem. 100, 1047-1056 (2004)
    J. Chem. Phys. 120, 8425 (2004)
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

(* erf_gau
    Int. J. of Quant. Chem. 100, 1047-1056 (2004)
*)
att_gau_aux1 := a -> sqrt(Pi)*erf(sqrt(3)/(2*a)):
att_gau_aux2 := a -> exp(-3/(4*a^2)):
attenuation_erf_gau := a ->
  8/3*a*(att_gau_aux1(a) - 2*sqrt(3)*a + 16*a^3/sqrt(27) + (2*a/sqrt(3) - 16*a^3/sqrt(27))*att_gau_aux2(a)):

(* erf_gau2
    J. Chem. Phys. 127, 154109 (2007)
    You can recover attenuation_erf_gau by putting a1 = 3
*)
att_gau2_aux1 := (a, a1) -> sqrt(Pi)*erf(sqrt(a1)/(2*a)):
att_gau2_aux2 := (a, a1) -> exp(-a1/(4*a^2)) - 1:
attenuation_gau2 := (a, a1) ->
  8/3*a*(att_gau2_aux1(a, a1) + 2*a/sqrt(a1)*att_gau2_aux2(a, a1)*(1 - 8*a^2/a1) - 4*a/sqrt(a1)):

(* yukawa
    Chem. Phys. Lett. 462(2008) 348-351
*)
att_yuk_aux1 := a -> arctan(1, a):
att_yuk_aux2 := a -> log(1.0 + 1.0/a^2):
att_yuk_aux3 := a -> a^2 + 1:
attenuation_yukawa := a -> my_piecewise3(
    a > 50, 1/(9*a^2) - 1/(30*a^4),
    1.0 - 8/3*a*(att_yuk_aux1(a) + a/4*(1 - (att_yuk_aux3(a) + 2)*att_yuk_aux2(a)))
):

