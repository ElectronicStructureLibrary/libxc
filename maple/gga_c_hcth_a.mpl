(* type: work_gga_c *)

gamma_ss := 0.2:
cc_ss    := [0.0136823, 0.268920, -0.550769,  1.03947, 0.0]:

gamma_ab := 0.006:
cc_ab    := [0.836897,  1.72051,  -2.78498,  -4.57504, 0.0]:

$include "lda_c_vwn.mpl"
$include "b97.mpl"

f := (rs, z, xt, xs0, xs1) ->
  f_b97(f_vwn, gamma_ss, cc_ss, gamma_ab, cc_ab,
        rs, z, xs0, xs1):