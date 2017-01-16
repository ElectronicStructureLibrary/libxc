(* type: work_mgga_x *)

gamm := (x, z) -> 1.0 + alpha*(x^2 + z):

gtv4 := (x, z) ->                                 \
  coeff_d[1]/gamm(x, z) +                               \
  (coeff_d[2]*x^2 + coeff_d[3]*z)/gamm(x,z)^2 +               \
  (coeff_d[4]*x^4 + coeff_d[5]*x^2*z + coeff_d[6]*z^2)/gamm(x,z)^3:
