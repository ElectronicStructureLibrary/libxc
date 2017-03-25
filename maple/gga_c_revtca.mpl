(* type: work_gga_c *)

$include "gga_c_tca.mpl"

msinc := x -> piecewise(x = 0, 1, sin(x)/x):
aa    := Pi*(9*Pi/4)^(1/3):

fD := (rs, z, s) -> 1 - z^4*(1 - msinc(aa*s/rs)^2):

f := (rs, z, xt, xs0, xs1) -> 
  f_tcs(rs, z, xt)*fD(rs, z, X2S*2^(1/3)*xt):