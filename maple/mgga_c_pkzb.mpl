(* type: work_mgga_c *)

$define gga_c_pbe_params
$include "gga_c_pbe.mpl"

pkzb_c := 0.53:

pkzb_perp := (rs, z, xt, xs0, xs1, ts0, ts1) ->
  + (1.0 + pkzb_c*(t_total(z, xs0^2, xs1^2)/(8.0*t_total(z, ts0, ts1)))^2)
  * f_pbe(rs, z, xt, xs0, xs1):

pkzb_par  := (rs, z, xt, xs0, xs1, ts0, ts1) ->
  - (1.0 + pkzb_c)*(
    + (xs0^2/(8.0*ts0))^2*gga_stoll_par(f_pbe, rs,  z, xs0,  1.0)
    + (xs1^2/(8.0*ts1))^2*gga_stoll_par(f_pbe, rs, -z, xs1, -1.0)
  ):

f_pkzb := (rs, z, xt, xs0, xs1, ts0, ts1) ->
  pkzb_perp(rs, z, xt, xs0, xs1, ts0, ts1) + pkzb_par(rs, z, xt, xs0, xs1, ts0, ts1):

f := (rs, z, xt, xs0, xs1, ts0, ts1, us0, us1) ->
  f_pkzb(rs, z, xt, xs0, xs1, ts0, ts1):
