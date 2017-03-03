(* type: work_gga_c *)

a := -0.74860:
b :=  0.06001:
c :=  3.60073:
d :=  0.90000:

f_num := (z, xt) -> sqrt(1.0 - z^2)*(a + b*xt):
f_den := (rs, xs0, xs1) -> c + d*(xs0 + xs1) + rs:

f := (rs, zeta, xt, xs0, xs1) -> f_num(zeta, xt)/f_den(rs, xs0, xs1):
