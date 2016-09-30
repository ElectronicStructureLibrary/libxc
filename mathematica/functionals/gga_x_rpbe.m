(* type: work_gga_x *)
(* prefix:
  gga_x_rpbe_params *params;
 
  assert(p->params != NULL);
  params = (gga_x_rpbe_params * )(p->params);
*)

f0[x_] = 1.0 + paramskappa (1.0 - Exp[-paramsmu x^2/paramskappa])
f[x_]  = f0[X2S x]
