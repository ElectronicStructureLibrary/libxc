(* type: work_gga_x *)

mu    := 0.804:
kappa := 0.219:
m     := 100:

gamm  := m*mu/kappa:
Cx    := kappa/m:

f0 := s -> 1.0 + sum(Cx * (gamm*s^2/(1.0 + gamm*s^2))^i, i=1..m):

f  := x -> f0(X2S*x):