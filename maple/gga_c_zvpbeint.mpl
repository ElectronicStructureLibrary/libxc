(* type: work_gga_c *)
(* prefix:
  gga_c_zvpbeint_params *params;

  assert(p->params != NULL);
  params = (gga_c_zvpbeint_params * )(p->params);
*)

params_a_gamma := (1.0 - log(2.0))/Pi^2:
params_a_BB    := 1.0:
$include "gga_c_pbe.mpl"

nu := (rs, z, t) ->
  t*mphi(z)*(3.0/rs)^(1/6):
ff := (rs, z, t) ->
  exp(-params_a_alpha*nu(rs, z, t)^3*m_abs(1.0*z)^params_a_omega):

f  := (rs, z, xt, xs0, xs1) ->
  f_pw(rs, z) + ff(rs, z, tp(rs, z, xt))*fH(rs, z, tp(rs, z, xt)):
