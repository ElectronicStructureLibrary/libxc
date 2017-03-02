(* type: work_mgga_c *)
(* prefix:
  mgga_c_bc95_params *params;

  assert(p->params != NULL);
  params = (mgga_c_bc95_params * )(p->params);
*)

$define lda_c_pw_params
$define lda_c_pw_modified_params
$include "lda_c_pw.mpl"

(* The B97 function g *)
bc95_gpar  := (xs, ts) -> (ts - xs^2/8.0)/(K_FACTOR_C*(1.0 + params_a_css*xs^2)^2):
bc95_gperp := (xs_0_, xs_1_) -> 1.0/(1.0 + params_a_copp*(xs_0_^2 + xs_1_^2)):

(* The parallel and perpendicular components of the energy *)
bc95_fpar  := (rs, z, xs_0_, xs_1_, ts_0_, ts_1_) ->
  + lda_stoll_par(f_pw, rs,  z,  1.0) * bc95_gpar(xs_0_, ts_0_)
  + lda_stoll_par(f_pw, rs, -z, -1.0) * bc95_gpar(xs_1_, ts_1_):

bc95_fperp := (rs, z, xs_0_, xs_1_) ->
  lda_stoll_perp(f_pw, rs, z) * bc95_gperp(xs_0_, xs_1_):

f_bc95 := (rs, z, xs_0_, xs_1_, ts_0_, ts_1_) ->
  + bc95_fpar (rs, z, xs_0_, xs_1_, ts_0_, ts_1_)
  + bc95_fperp(rs, z, xs_0_, xs_1_):

f := (rs, z, xt, xs_0_, xs_1_, ts_0_, ts_1_, us_0_, us_1_) ->
 f_bc95(rs, z, xs_0_, xs_1_, ts_0_, ts_1_):

