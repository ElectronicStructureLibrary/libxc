(* type: work_mgga_c *)

zlp_c := 0.828432   *RS_FACTOR:
zlp_d := 2.15509e-2 *RS_FACTOR:
zlp_k := 2.047107e-3*RS_FACTOR:

f_zlp := (rs, z, xt, us0, us1) ->
  - (zlp_c + zlp_d*t_vw(z, xt, us0, us1))
  * (1 - zlp_k*log(1 + rs/zlp_k)/rs)/rs:

f := (rs, z, xt, xs0, xs1, ts0, ts1, us0, us1) ->
  f_zlp(rs, z, xt, us0, us1):
