(* type: work_lda *)

C1 := -0.0603:
C2 :=  0.0175:
C3 := -0.00053:

f := (rs, zeta) -> C1 + C2*n_total(rs)^(-1/3) + C3*n_total(rs)^(-2/3):

