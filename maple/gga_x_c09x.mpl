(* type: work_gga_x *)

mu    := 0.0617:
kappa := 1.245:
alpha := 0.0483:

f0 := s -> 1.0 + mu*s^2*exp(-alpha*s^2) + kappa*(1.0 - exp(-0.5*alpha*s^2)):
f  := x -> f0(X2S*x):
