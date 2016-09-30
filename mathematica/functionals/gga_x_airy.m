(* type: work_gga_x *)

f1[s_] = (1 - a5 s^a6 + a7 s^a8)/(1 + a9 s^a10)

Import["gga_x_lag.m"]
Clear[f]

f[x_] = f0[X2S*x] + f1[X2S*x]
