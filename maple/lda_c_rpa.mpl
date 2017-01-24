(* type: work_lda *)

a :=  0.0311:
b := -0.048:
c :=  0.009:
d := -0.017:

f := (rs, zeta) -> a*log(rs) + b + c*rs*log(rs) + d*rs:
