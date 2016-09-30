(* type: work_gga_x *)

f0[s_] = 1.0 + B1 s Log[1.0 + s] + B2 s Log[1.0 + Log[1.0 + s]]

f[x_] = f0[X2S*x]
