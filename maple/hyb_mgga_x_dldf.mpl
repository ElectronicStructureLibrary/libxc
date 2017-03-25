(* type: work_mgga_x *)

params_a_kappa := 4.8827323:
params_a_mu    := 0.3511128:
$include "gga_x_pbe.mpl"

dldf_a := [1, -0.1637571, -0.1880028, -0.4490609, -0.0082359]:
csi_HF := 1 - 0.6144129:

f := (rs, x, t, u) ->
  + csi_HF*f_pbe(x)*mgga_series_w(dldf_a, 5, t):
