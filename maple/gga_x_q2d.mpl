(* type: work_gga_x *)

(* standard PBE_SOL *)
kappa := 0.8040:
mu    := MU_GE:
$include "pbe_x.mpl"

cc := 100.0:
c1 := 0.5217:

f1 := s -> pbe_x(s/X2S)*(cc - s^4) + c1*s^3.5*(1.0 + s^2):
f2 := s -> cc + s^6:

f  := x -> f1(X2S*x)/f2(X2S*x):