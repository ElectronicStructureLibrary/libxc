(* type: work_gga_x *)

f0 := s -> 1.0 + B1*s*log(1.0 + s) + B2*s*log(1.0 + log(1.0 + s)):
f  := x -> f0(X2S*x):
