(* type: work_lda *)

$include "vwn.mpl"

f := (rs, z) ->
  + f_aux(A_vwn[1], b_vwn[1], c_vwn[1], x0_vwn[1], rs)
  + f_aux(A_rpa[3], b_rpa[3], c_rpa[3], x0_rpa[3], rs)*f_zeta(z)*(1 - z^4)/fpp
  - DRPA(rs, z)*f_zeta(z)*(1 - z^4)
  +  DMC(rs, z)*f_zeta(z):