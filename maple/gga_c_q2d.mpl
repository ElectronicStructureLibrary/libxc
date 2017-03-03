(* type: work_gga_c *)

$define gga_c_pbe_params
$include "gga_c_pbe.mpl"

$include "lda_c_2d_amgb.mpl"

rs2D_factor := 1.704:
q2d_dd := 1e6:

q2d_rs2D := (rs, xt) -> rs2D_factor*rs*sqrt(X2S*xt)/RS_FACTOR:

fac := t -> t^4*(1.0 + t^2)/(q2d_dd + t^6):

f_q2d := (rs, z, xt, xs0, xs1) ->
  (1.0 - fac(tt(rs, z, xt)))*f_pbe(rs, z, xt, xs0, xs1) + fac(tt(rs, z, xt))*f_amgb(q2d_rs2D(rs, xt), z):

f := (rs, z, xt, xs0, xs1) -> f_q2d(rs, z, xt, xs0, xs1):
