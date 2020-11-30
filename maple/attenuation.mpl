(*
 Copyright (C) 2017 M.A.L. Marques
               2020 Susi Lehtola

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
(* This is the full function, which is numerically unstable for large a *)
attenuation_erf0 := a ->
  1 - 8/3*a*(att_erf_aux1(a) + 2*a*(att_erf_aux2(a) - att_erf_aux3(a))):
(* For small a, we use a Taylor expansion at 1/a=0 since that's the
   problem in aux1 and aux2.  For some reason, Maple 2020 doesn't
   compute the Taylor series up to 6th order in a unless one specifies
   12th order for the evaluation...  *)
   attenuation_erf_taylor := a ->
eval(convert(taylor(eval(attenuation_erf0(b), b=1/inva), inva=0, 12), polynom), inva=1/a):
(* Now we just glue these together. The Taylor series is roughly
        1/36 * inva^2 - 1/960 * inva^4 + 1/26880 inva^6 + O(inva^8)
      = 1/36 * inva^2 * ( 1 - 3/80 * inva^2 + 3/2240 * inva^4 + O(inva^6))
so we can safely glue the pieces together at a = DBL_EPSILON^(-1/4)
since the truncated Taylor is numerically exact *)
attenuation_erf := a -> my_piecewise3( a >= DBL_EPSILON^(-1/4), attenuation_erf_taylor(a), attenuation_erf0(a)):

(* These are for hyb_mgga_x_js18 and hyb_mgga_x_pjs18. The terms that go to zero for large a: *)
attenuation_erf_f2_aux0 := a -> (20*a^2 - 64*a^4)*exp(-1/(4*a^2)) + 10*sqrt(Pi)*a*erf(1/(2*a)):
(* To make this go to zero for large a, we need to throw out the terms with negative powers of 1/a from the series expansion *)
throw_out_negative_inva := proc(X) select(t -> degree(t, {inva})>0, X); end proc:
attenuation_erf_f2_aux_taylor := a ->
eval(throw_out_negative_inva(convert(series(eval(attenuation_erf_f2_aux0(b), b=1/inva), inva=0, 12), polynom)), inva=1/a):
attenuation_erf_f2_aux := a -> my_piecewise3( a >= DBL_EPSILON^(-1/4), attenuation_erf_f2_aux_taylor(a), attenuation_erf_f2_aux0(a)):
attenuation_erf_f2 := a ->
  1 + 24*a^2*(- 3 - 36*a^2 + 64*a^4 + attenuation_erf_f2_aux(a)):

attenuation_erf_f3_aux0 := a -> (-8*a + 256*a^3 - 576*a^5 + 3849*a^7 - 122880*a^9)*exp(-1/(4*a^2)) + 2*sqrt(Pi)*(-2 + 60*a^2)*erf(1/(2*a)):
attenuation_erf_f3_aux_taylor := a ->
eval(throw_out_negative_inva(convert(series(eval(attenuation_erf_f3_aux0(b), b=1/inva), inva=0, 12), polynom)), inva=1/a):
attenuation_erf_f3_aux := a -> my_piecewise3( a >= DBL_EPSILON^(-1/4), attenuation_erf_f3_aux_taylor(a), attenuation_erf_f3_aux0(a)):
attenuation_erf_f3 := a ->
  1 + 8*a/7*(
      + 24*a^3*(-35 + 224*a^2 - 1440*a^4 + 5120*a^6)
      + attenuation_erf_f3_aux(a)):

(* erf_gau - screening function = + 2 mu/sqrt(pi) exp(-mu^2 r^2)
    Song et al, J. Chem. Phys. 127, 154109 (2007); doi:10.1063/1.2790017
    You can recover the result in Int. J. of Quant. Chem. 100, 1047 (2004)
    by putting a = a/sqrt(3) and multiplying the whole attenuation function by -sqrt(3)
*)
attenuation_gau := (a) ->
  -8/3*a*(att_erf_aux1(a) + 2*a*att_erf_aux2(a)*(1 - 8*a^2) - 4*a):

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

