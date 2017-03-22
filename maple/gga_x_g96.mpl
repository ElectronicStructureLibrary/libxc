(* type: work_gga_x *)

g96_c1 := 1.0/137.0:

f_g96  := x->
  1.0 + g96_c1/X_FACTOR_C * x^(3.0/2.0):

f  := x -> f_g96(x):
