(* type: work_lda *)

f_zeta_k := z -> 1/2*((1 + z)^(5/3) + (1 - z)^(5/3)):

c1 := 3.2372*RS_FACTOR:
c2 := 0.00196*RS_FACTOR:

f := (rs, zeta) -> c1*f_zeta_k(zeta)/rs^2
  * (1 - c2/rs*log(1 + rs/c2)):