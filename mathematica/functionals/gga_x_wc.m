(* type: work_gga_x *)

f0[s_]  = kappa + MUGE s^2 + (mu - MUGE) s^2 Exp[-s^2] + Log[1.0 + c s^4]

f[x_] = 1.0 + kappa*(1.0 - kappa/f0[X2S*x])
