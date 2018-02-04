(* type: work_mgga_c *)

lp90_c0 := 0.80569:
lp90_d0 := 3.0124e-3:
lp90_k  := 4.0743e-3:

(* Equation (60) *)
f_lp90 := (rs, z, xt, us0, us1) ->
  - (lp90_c0 + lp90_d0*t_vw(z, xt, us0, us1))/(rs/RS_FACTOR + lp90_k):

f := (rs, z, xt, xs0, xs1, ts0, ts1, us0, us1) ->
  f_lp90(rs, z, xt, us0, us1):
