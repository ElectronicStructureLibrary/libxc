(* type: work_gga_c *)

$include "gga_c_tca.mpl"

msinc := x -> piecewise(x = 0, 1.0, sin(x)/x):
aa    := Pi*(9.0*Pi/4.0)^(1.0/3.0):

fD := (rs, z, s) -> 1.0 - z^4*(1.0 - msinc(aa*s/rs)^2):

f := (rs, z, xt, xs_0_, xs_1_) -> 
  f_tcs(rs, z, xt)*fD(rs, z, X2S*2.0^(1.0/3.0)*xt):