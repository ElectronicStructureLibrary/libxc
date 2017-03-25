(* type: work_lda *)

a1 := -0.0357:
a2 :=  0.0562:
b1 := -0.0311:
b2 :=  2.39:

f := (rs, zeta) -> a1/(1 + a2*rs/RS_FACTOR)
  + b1*log((rs/RS_FACTOR + b2)/(rs/RS_FACTOR)):
