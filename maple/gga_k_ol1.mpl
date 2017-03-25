(* type: work_gga_x *)

c4 := 0.00677:

f0 := s -> 1 + (s^2/72 + c4*s)/K_FACTOR_C:
f  := x -> f0(x/2^(1/3)):