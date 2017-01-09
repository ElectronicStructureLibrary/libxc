(* type: work_gga_x *)

f0 := s -> kappa + MU_GE * s^2 + (mu - MU_GE) * s^2 * exp(-s^2) + log(1.0 + c * s^4):

f  := x -> 1.0 + kappa*(1.0 - kappa/f0(X2S*x)):
