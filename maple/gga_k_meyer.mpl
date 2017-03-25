(* type: work_gga_x *)

lambda := y -> 1/2*(1 + (1 - y^2)*log((1 + y)/abs(1 - y))/(2*y)):

f := x -> 1 + lambda(X2S*x/6)*x^2/(8*K_FACTOR_C):