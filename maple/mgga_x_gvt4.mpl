(* type: work_mgga_x *)

alpha := 0.00186726:
d     := [-9.800683e-01, -3.556788e-03, 6.250326e-03, -2.354518e-05, -1.282732e-04, 3.574822e-04]:

gamm := (x, z) -> 1.0 + alpha*(x^2 + z):

gtv4 := (x, z) ->                                \
  d[1]/gamm(x, z) +                               \
  (d[2]*x^2 + d[3]*z)/gamm(x,z)^2 +               \
  (d[4]*x^4 + d[5]*x^2*z + d[6]*z^2)/gamm(x,z)^3:

f := (rs, x, t, u) -> -gtv4(x, 2.0*(t - K_FACTOR_C))/X_FACTOR_C: