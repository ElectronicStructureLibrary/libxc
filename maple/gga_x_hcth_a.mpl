(* type: work_gga_x *)

mbeta  :=  0.0042:
mgamma :=  6:
c0    :=  1.09878:
c1    := -2.51173:
c2    :=  0.0156233:

f_aux := x -> 1.0 + mgamma*mbeta*x*arcsinh(x):

f := x -> c0 + mbeta/X_FACTOR_C*x^2*(c1/f_aux(x) + c2/(mbeta*f_aux(x)^2)):