(* type: work_gga_c *)

$include "lda_c_rc04.mpl"

msigma := 1.43:
malpha := 2.30:

Bs := s -> 1.0/(1.0 + msigma*s^malpha):

f_tcs := (rs, z, xt) -> f_rc04(rs, z)*Bs(X2S*2.0^(1.0/3.0)*xt):

f := (rs, z, xt, xs_0_, xs_1_) -> 
  f_tcs(rs, z, xt):