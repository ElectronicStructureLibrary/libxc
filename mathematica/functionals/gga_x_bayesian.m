(* type: work_gga_x *)

f0[s_] = s^2/(1.0 + s)^2

f[x_] = theta0 + f0[X2S*x] (theta1 + f0[X2S*x] theta2)
