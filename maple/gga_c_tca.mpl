(* type: work_gga_c *)

$include "lda_c_rc04.mpl"

msigma := 1.43:
malpha := 2.30:

Bs := s -> 1/(1 + msigma*s^malpha):

f_tcs := (rs, z, xt) -> f_rc04(rs, z)*Bs(X2S*2^(1/3)*xt):

f := (rs, z, xt, xs0, xs1) -> 
  f_tcs(rs, z, xt):