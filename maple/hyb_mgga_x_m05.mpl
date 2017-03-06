(* type: work_mgga_x *)
(* prefix:
  mgga_x_m05_params *params;

  assert(pt->params != NULL);
  params = (mgga_x_m05_params * )(pt->params);
*)

$define gga_x_pbe_params
$include "gga_x_pbe.mpl"

f := (rs, x, t, u) ->
  + params_a_csi_HF*f_pbe(x)*mgga_series_w(params_a_a, 12, t):
