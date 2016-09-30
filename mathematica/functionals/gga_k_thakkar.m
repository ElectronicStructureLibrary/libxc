(* type: work_gga_x *)
(* prefix:
  gga_x_b88_params *params;
 
  assert(p->params != NULL);
  params = (gga_x_b88_params * )(p->params);
*)

Import["gga_x_b88.m"]
Clear[f]

f1[x_] = -0.072*x/(1.0 + 2.0*4.0^(1.0/3.0) x)

f[x_] = f0[x] + f1[x]
