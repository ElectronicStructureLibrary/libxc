(* type: work_lda *)

f_zeta_k := z -> 0.5*((1.0 + z)^(5.0/3.0) + (1.0 - z)^(5.0/3.0)):

c1 := 3.2372*RS_FACTOR:
c2 := 0.00196*RS_FACTOR:

f := (rs, zeta) -> c1*f_zeta_k(zeta)/rs^2
  * (1.0 - c2/rs*log(1.0 + rs/c2)):