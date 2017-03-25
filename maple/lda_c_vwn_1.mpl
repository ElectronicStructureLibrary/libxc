(* type: work_lda *)

$include "vwn.mpl"

f := (rs, z) -> 
  + f_aux(A_vwn[1], b_vwn[1], c_vwn[1], x0_vwn[1], rs)*(1 - f_zeta(z)) 
  + f_aux(A_vwn[2], b_vwn[2], c_vwn[2], x0_vwn[2], rs)*f_zeta(z):