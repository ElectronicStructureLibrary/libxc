(* type: work_gga_x *)
(* prefix:
  gga_x_sogga11_params *params;
 
  assert(p->params != NULL);
  params = (gga_x_sogga11_params * )(p->params);
*)

alpha := params_a_mu*X2S*X2S/params_a_kappa:

f0 := x -> 1.0 - 1.0/(1.0 + alpha*x^2):
f1 := x -> 1.0 - exp(-alpha*x^2):

f  := x -> add(params_a_a[i]*f0(x)^(i-1), i=1..6) + add(params_a_b[i]*f1(x)^(i-1), i=1..6):