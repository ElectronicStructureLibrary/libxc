(* type: work_lda *)

a0 := 0.93222*RS_FACTOR:
kk := 9.47362e-3*RS_FACTOR:

f := (rs, zeta) -> -a0*(1.0 - kk*log(1.0 + rs/kk)/rs)/rs:
