(* type: work_mgga_x *)
(* prefix:
  mgga_x_tau_hcth_params *params;

  assert(pt->params != NULL);
  params = (mgga_x_tau_hcth_params * ) (pt->params);
*)

coeff_a := [0, 1, 0, -2, 0, 1]:

(* Equation (29) *)
gamX := 0.004:
ux   := x -> gamX*x^2/(1 + gamX*x^2):

gxl  := x -> add(params_a_cx_local [i]*ux(x)^(i-1), i=1..4):
gxnl := x -> add(params_a_cx_nlocal[i]*ux(x)^(i-1), i=1..4):

f    := (rs, x, t, u) -> gxl(x) + gxnl(x)*mgga_series_w(coeff_a, 6, t):
