(* type: work_gga_x *)

kappa := 0.8040:
mu    := 0.2195149727645171:
m     := 100:

gamm  := m*mu/kappa:
Cx    := kappa/m:

f0 := s -> 1 + add(Cx * (gamm*s^2/(1 + gamm*s^2))^i, i=1..m):

f  := x -> f0(X2S*x):
