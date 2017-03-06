(* type: work_mgga_c *)
(* prefix:
  mgga_c_m06l_params *params;

  assert(p->params != NULL);
  params = (mgga_c_m06l_params * )(p->params);
*)

$include "mgga_c_vsxc.mpl"
$include "mgga_c_m05.mpl"

f_m06 := (rs, z, xs0, xs1, ts0, ts1) ->
  + f_m05(rs, z, xs0, xs1, ts0, ts1)
  + f_vsxc(rs, z, xs0, xs1, ts0, ts1):

f := (rs, z, xt, xs0, xs1, ts0, ts1, us0, us1) ->
  f_m06(rs, z, xs0, xs1, ts0, ts1):

