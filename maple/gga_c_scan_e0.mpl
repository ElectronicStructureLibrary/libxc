(* type: work_gga_c *)

$define gga_c_pbe_params
$include "gga_c_regtpss.mpl"

f1 := (rs, z, t) -> 1 + 4*A(rs, z, t)*t^2:
f2 := (rs, z, t) -> mbeta(rs, t)*(1 - f1(rs, z, t)^(-1/4))/(mgamma*A(rs, z, t)):