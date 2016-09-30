(* type: work_gga_x *)

csi[s_] = (3.0/2.0 ProductLog[s^(3.0/2.0) / (2.0 Sqrt[6.0])])^(2.0/3.0)

fb[s_]  = Pi/3.0 s/(csi[s] (d + csi[s]^2)^(1.0/4.0))

flaa[s_]  = (1.0 + c s^2)/(1.0 + c s^2/fb[s])

XX[s_] = 1.0 - alpha s^2/(1.0 + alpha s^2)

f[x_] = XX[X2S*x] + (1.0 - XX[X2S*x])*flaa[X2S*x]
