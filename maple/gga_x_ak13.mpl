(* type: work_gga_x *)

f0 := s -> 1 + B1*s*log(1 + s) + B2*s*log(1 + log(1 + s)):
f  := x -> f0(X2S*x):
