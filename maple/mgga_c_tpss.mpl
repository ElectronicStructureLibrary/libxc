(* type: work_mgga_c *)
(* prefix:
  mgga_c_tpss_params *params;

  assert(p->params != NULL);
  params = (mgga_c_tpss_params * )(p->params);
*)

(* beta is taken from the params *)
params_a_gamma := (1.0 - log(2.0))/Pi^2:
params_a_BB    := 1.0:
$include "gga_c_pbe.mpl"

$include "tpss.mpl"

f := (rs, z, xt, xs0, xs1, ts0, ts1, us0, us1) ->
  + f_tpss(f_pbe, rs, z, xt, xs0, xs1, ts0, ts1):




