(* type: work_mgga_c *)

$include "gga_c_regtpss.mpl"

params_a_C0_c := [0.59, 0.9269, 0.6225, 2.1540]:
params_a_d    := 2.8:
$include "tpss.mpl"

f := (rs, z, xt, xs0, xs1, ts0, ts1, us0, us1) ->
  + f_tpss(f_pbe, rs, z, xt, xs0, xs1, ts0, ts1):




