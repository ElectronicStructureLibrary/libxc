(* type: work_gga_x *)

$define gga_x_pbe_tca_params
$include "gga_x_pbe.mpl"

$define gga_x_pw91_params
$include "gga_x_pw91.mpl"

malpha :=  1.0:
mbeta  := 19.0:

fab := x -> 1.0/(1.0 + exp(-malpha*(x - mbeta))):
f   := x -> (1.0 - fab(x))*f_pbe(x) + fab(x)*f_pw91(x):