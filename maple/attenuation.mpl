(*
 Copyright (C) 2017 M.A.L. Marques
               2020 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)


(* Cap polynomial expansion to given order *)
throw_out_large_n := proc(X,n) select(t -> abs(degree(t,{b}))<=n, X); end proc:

(* Function that makes f(a) smooth for a->infty *)
enforce_smooth_lr := proc(f, a);

  (* The order to use for the expansion: 80th order.

     I know this sounds ridiculously high, but it turns out to be
     necessary to avoid numerical difficulties in the *original*
     function.

     We want to make the original function f(a) and its series
     expansion f_series(a) match numerically exactly at the cutoff
     point, but numerical experiments (and failing test suites)
     demonstrate f(a) is unstable unless a very high-order expansion
     is used, since it allows bringing the cutoff point to lower
     values of a where the functions are more stable.

     Much lower-order expansions could be used if we had e.g. 64-digit
     accuracy (instead of the ~16 digits with double precision), but
     this is what we have to live with. Note, however, that since the
     expansions appear to only use (1/a)^2, the actual number of terms
     in the Taylor series is half of this...

     erf and Yukawa are fine at order 50, but f30 is still too
     unstable and needs 80...

     2020-12-07 Susi Lehtola
  *)
  local expansion_order := 80:

  (* Due to a bug in Maple 2020, series aren't being computed to the
     requested order. So, we need to use a larger expansion order and
     then truncate it back. Pad the expansion by this order *)
   local padding_order := 30:

  (* Calculate large-a expansion *)
  f_large := a -> eval(throw_out_large_n(convert(series(f(b), b=infinity, expansion_order+padding_order), polynom), expansion_order), b=a):

  (* The first term is usually inva^2, which gives an expansion of the type
             inva^2(O(1) + O(inva^2) + O(inva^4) + ...)
     meaning that the last term in the bracket is O(inva^{expansion_order-2}).
     This term is numerically zero for a larger than
  *)
  local a_cutoff := DBL_EPSILON^(-1/(expansion_order-2)):
  (* Return the series expansion for large a; also remove any numerical overflows from the original branch  *)
  my_piecewise3(a >= a_cutoff, f_large(m_max(a, a_cutoff)), f(m_min(a, a_cutoff))):
end proc:

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
attenuation_erf := a -> enforce_smooth_lr(attenuation_erf0, a):

(* These are for hyb_mgga_x_js18 and hyb_mgga_x_pjs18. This is the
bracket in eqn (10) in Patra et al, 2018 *)
attenuation_erf_f20 := a ->
  1 + 24*a^2*( (20*a^2 - 64*a^4)*exp(-1/(4*a^2))
    - 3 - 36*a^2 + 64*a^4
    + 10*a*sqrt(Pi)*erf(1/(2*a))):
attenuation_erf_f2 := a -> enforce_smooth_lr(attenuation_erf_f20, a):

(* This is eqn (11) in Patra et al, 2018  *)
attenuation_erf_f30 := a ->
  1 + 8/7*a*(
        (-8*a + 256*a^3 - 576*a^5 + 3840*a^7 - 122880*a^9)*exp(-1/(4*a^2))
      + 24*a^3*(-35 + 224*a^2 - 1440*a^4 + 5120*a^6)
      + 2*sqrt(Pi)*(-2 + 60*a^2)*erf(1/(2*a))):
attenuation_erf_f3 := a -> enforce_smooth_lr(attenuation_erf_f30, a):

(* erf_gau - screening function = + 2 mu/sqrt(pi) exp(-mu^2 r^2)
    Song et al, J. Chem. Phys. 127, 154109 (2007); doi:10.1063/1.2790017
    You can recover the result in Int. J. of Quant. Chem. 100, 1047 (2004)
    by putting a = a/sqrt(3) and multiplying the whole attenuation function by -sqrt(3)
*)
attenuation_gau0 := a -> -8/3*a*(att_erf_aux1(a) + 2*a*att_erf_aux2(a)*(1 - 8*a^2) - 4*a):
attenuation_gau := a -> enforce_smooth_lr(attenuation_gau0, a):

(* yukawa
    Akinaga and Ten-no, Chem. Phys. Lett. 462, 348 (2008); doi:10.1016/j.cplett.2008.07.103
*)
att_yuk_aux1 := a -> arctan(1, a):
att_yuk_aux2 := a -> log(1 + 1/a^2):
att_yuk_aux3 := a -> a^2 + 1:
attenuation_yukawa0 := a -> 1 - 8/3*a*(att_yuk_aux1(a) + a/4*(1 - (att_yuk_aux3(a) + 2)*att_yuk_aux2(a))):
attenuation_yukawa := a -> enforce_smooth_lr(attenuation_yukawa0, a):
