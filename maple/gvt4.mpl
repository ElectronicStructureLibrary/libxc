(* type: work_mgga_x *)

gvt4_gamm := (alpha, x, z) -> 1 + alpha*(x^2 + z):

gtv4 := (alpha, dd, x, z) ->
  dd[1]/gvt4_gamm(alpha, x, z) +
  (dd[2]*x^2 + dd[3]*z)/gvt4_gamm(alpha, x, z)^2 +
  (dd[4]*x^4 + dd[5]*x^2*z + dd[6]*z^2)/gvt4_gamm(alpha, x, z)^3:
