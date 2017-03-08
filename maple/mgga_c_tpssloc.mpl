(* type: work_mgga_c *)

$include "gga_c_pbeloc.mpl"

params_a_C0_c := [0.35, 0.87, 0.50, 2.26]:
params_a_d    := 4.5:
$include "tpss.mpl"

f := (rs, z, xt, xs0, xs1, ts0, ts1, us0, us1) ->
  + f_tpss(f_pbe, rs, z, xt, xs0, xs1, ts0, ts1):
