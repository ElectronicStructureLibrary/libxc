(* type: work_mgga_c *)

$define gga_x_b88_params
$include "gga_x_b88.mpl"

cab := 0.63:
css := 0.96:

b88_cab := (rs, z, xs0, xs1) ->
  - 0.8 * (1.0 - z^2)/4.0 * n_total(rs)
  * b88_zab(cab, f_b88, rs, z, xs0, xs1)^2 * (
    1.0 - log(1.0 + b88_zab(cab, f_b88, rs, z, xs0, xs1))
      / b88_zab(cab, f_b88, rs, z, xs0, xs1)
    ):

b88_css := (rs, z, xs, ts) ->
  - 0.01 * (1.0 + z)/2.0 * n_spin(rs, z)^(5/3) * (2.0*ts - xs^2/4.0)
  * b88_zss(css, f_b88, rs, z, xs)^4 * (
    1.0 - 2.0*log(1.0 + b88_zss(css, f_b88, rs, z, xs)/2.0)
      / b88_zss(css, f_b88, rs, z, xs)
    ):

f_b88_c := (rs, z, xs0, xs1, ts0, ts1) ->
  + b88_cab(rs,  z, xs0, xs1)
  + b88_css(rs,  z, xs0, ts0)
  + b88_css(rs, -z, xs1, ts1):

f := (rs, z, xt, xs0, xs1, ts0, ts1, us0, us1) ->
  f_b88_c(rs, z, xs0, xs1, ts0, ts1):
