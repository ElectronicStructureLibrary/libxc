(* type: work_mgga_c *)

$include "gga_c_gapc.mpl"

(* override definition of gap_C *)
gap_G := (rs, z, xt, par) -> RS_FACTOR^2/8.0 * xt^2/rs^2:

f_kcis := (rs, z, xt, xs0, xs1, ts0, ts1) ->
  + f_gap(rs, z, xt)
  - xs0^2/(8.0*ts0) * f_gap(rs,  1.0, xs0)
  - xs1^2/(8.0*ts1) * f_gap(rs, -1.0, xs1):

f  := (rs, z, xt, xs0, xs1, ts0, ts1, us0, us1) ->
  f_kcis(rs, z, xt, xs0, xs1, ts0, ts1):
