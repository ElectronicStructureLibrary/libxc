(* type: work_gga_x *)

g96_c1 := 1/137:

f_g96  := x->
  1 + g96_c1/X_FACTOR_C * x^(3/2):

f  := x -> f_g96(x):
