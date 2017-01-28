(* type: work_gga_x *)

mbeta  := 0.003317:
mgamma := 0.008323:

f := x -> 1.0 + mbeta/X_FACTOR_C*x^2/(1.0 + mgamma*x^2)^(3.0/4.0):