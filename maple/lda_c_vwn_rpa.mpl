(* type: work_lda *)

$include "vwn.mpl"

f := (rs, z) ->
  + f_aux(A_rpa[1], b_rpa[1], c_rpa[1], x0_rpa[1], rs)*(1 - f_zeta(z)) 
  + f_aux(A_rpa[2], b_rpa[2], c_rpa[2], x0_rpa[2], rs)*f_zeta(z):