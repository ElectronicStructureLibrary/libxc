(* type: work_lda *)

AA := -0.655868:
BB :=  4.888270:
CC :=  3.177037:
DD :=  0.897889:

phi := z -> 0.5*((1.0 + z)^(2.0/3.0) + (1.0 - z)^(2.0/3.0)):

f_rc04 := (rs, zeta) -> phi(zeta)^3 * (AA*arctan(BB + CC*rs) + DD)/rs:

f      := (rs, zeta) -> f_rc04(rs, zeta):