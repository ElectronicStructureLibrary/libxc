(* type: work_gga_x *)

mbeta := 0.018641:
mcsi  := 8.0: (* for harmonic potentials *)

f := x -> 1.0 + mbeta/X_FACTOR_2D_C*x^2/(1.0 + mcsi*mbeta*x*arcsinh(x)):