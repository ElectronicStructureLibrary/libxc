(* type: work_mgga_c *)
(* prefix:
  mgga_c_m08_params *params;

  assert(p->params != NULL);
  params = (mgga_c_m08_params * )(p->params);
*)


$define lda_c_pw_params
$define lda_c_pw_modified_params
$include "lda_c_pw.mpl"

$define gga_c_pbe_params
$include "gga_c_pbe.mpl"

(* the prefactor was chosen to get the right K_FACTOR_C in mgga_series_w *)
f_m08 := (rs, z, xt, xs_0_, xs_1_, ts_0_, ts_1_) ->
  + mgga_series_w(params_a_m08_a, 12, 2.0^(2.0/3.0)*t_total(z, ts_0_, ts_1_))
    * f_pw(rs, z)
  + mgga_series_w(params_a_m08_b, 12, 2.0^(2.0/3.0)*t_total(z, ts_0_, ts_1_))
    * (f_pbe(rs, z, xt, xs_0_, xs_1_) - f_pw(rs, z)):

f := (rs, z, xt, xs_0_, xs_1_, ts_0_, ts_1_, us_0_, us_1_) ->
  f_m08(rs, z, xt, xs_0_, xs_1_, ts_0_, ts_1_, us_0_, us_1_):
