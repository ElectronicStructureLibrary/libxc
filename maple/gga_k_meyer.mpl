(* type: work_gga_x *)

lambda := y -> 0.5*(1.0 + (1.0 - y^2)*log((1.0 + y)/abs(1.0 - y))/(2.0*y)):

f := x -> 1.0 + lambda(X2S*x/6.0)*x^2/(8.0*K_FACTOR_C):