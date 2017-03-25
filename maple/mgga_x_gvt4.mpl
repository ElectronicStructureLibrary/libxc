(* type: work_mgga_x *)

malpha  := 0.00186726:
coeff_d := [-9.800683e-01, -3.556788e-03, 6.250326e-03, -2.354518e-05, -1.282732e-04, 3.574822e-04]:

$include "gvt4.mpl"

f := (rs, x, t, u) -> -gtv4(malpha, coeff_d, x, 2*(t - K_FACTOR_C))/X_FACTOR_C:
