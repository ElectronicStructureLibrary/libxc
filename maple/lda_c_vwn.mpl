(* type: work_lda *)

$include "vwn.mpl"

f_vwn := (rs, z) ->
  + f_aux(A_vwn[1], b_vwn[1], c_vwn[1], x0_vwn[1], rs)
  + f_aux(A_vwn[3], b_vwn[3], c_vwn[3], x0_vwn[3], rs)*f_zeta(z)*(1 - z^4)/fpp
  +  DMC(rs, z)*f_zeta(z)*z^4:

f := (rs, z) -> f_vwn(rs, z):