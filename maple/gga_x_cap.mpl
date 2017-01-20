(* type: work_gga_x *)

Ax       := -3.0/4.0*(3.0/Pi)^(1.0/3.0):
mu       := 0.2195149727645171:
alpha    := -Ax*mu:
c        := alpha/(3.0*Pi^2)^(1.0/3.0):
alphaoAx := -mu:

f0 := s -> 1.0 - alphaoAx*s*log(1.0 + s)/(1.0 + c*log(1.0 + s)):
f  := x -> f0(X2S*x):