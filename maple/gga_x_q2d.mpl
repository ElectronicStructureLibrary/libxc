(* type: work_gga_x *)

$define gga_x_pbe_sol_params
$include "gga_x_pbe.mpl"

cc := 100.0:
c1 := 0.5217:

f1 := s -> f0_pbe(s)*(cc - s^4) + c1*s^3.5*(1.0 + s^2):
f2 := s -> cc + s^6:

f  := x -> f1(X2S*x)/f2(X2S*x):