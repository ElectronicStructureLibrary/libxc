(* type: work_gga_c *)
(* prefix:
  gga_c_zpbeint_params *params;
 
  assert(p->params != NULL);
  params = (gga_c_zpbeint_params * )(p->params);
*)

params_a_gamma := (1.0 - log(2.0))/Pi^2:
params_a_BB    := 1.0:
$include "gga_c_pbe.mpl"

ff := (z, t) -> mphi(z)^(params_a_alpha*t^3):

f  := (rs, z, xt, xs0, xs1) ->
  f_pw(rs, z) + ff(z, tp(rs, z, xt))*fH(rs, z, tp(rs, z, xt)):