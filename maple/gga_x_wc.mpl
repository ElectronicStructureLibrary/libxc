(* type: work_gga_x *)

mu    := 0.2195149727645171:
c     := (146.0/2025.0)*(4.0/9.0) - (73.0/405.0)*(2.0/3.0) + (mu - 10.0/81.0):
kappa := 0.8040:

f0 := s -> kappa + MU_GE * s^2 + (mu - MU_GE) * s^2 * exp(-s^2) + log(1.0 + c * s^4):

f  := x -> 1.0 + kappa*(1.0 - kappa/f0(X2S*x)):