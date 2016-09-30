(* type: work_gga_x *)

f0[s_] = 1.0 + mu s^2 Exp[-alpha s^2] + kappa (1.0 - Exp[-0.5 alpha s^2])

f[x_] = f0[X2S x]
