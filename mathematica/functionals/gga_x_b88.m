(* type: work_gga_x *)
(* prefix:
  gga_x_b88_params *params;
 
  assert(p->params != NULL);
  params = (gga_x_b88_params * )(p->params);
*)

f0[x_] = 1. + paramsbeta/XFACTORC x^2/(1.0 + paramsgamma paramsbeta x ArcSinh[x])

f[x_] = f0[x]
