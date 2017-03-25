(* type: work_lda *)

AA := -0.655868:
BB :=  4.888270:
CC :=  3.177037:
DD :=  0.897889:

phi := z -> 1/2*((1 + z)^(2/3) + (1 - z)^(2/3)):

f_rc04 := (rs, zeta) -> phi(zeta)^3 * (AA*arctan(BB + CC*rs) + DD)/rs:

f      := (rs, zeta) -> f_rc04(rs, zeta):