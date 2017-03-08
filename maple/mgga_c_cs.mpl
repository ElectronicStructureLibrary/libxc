(* type: work_mgga_c *)

cs_a := -0.04918:
cs_b :=  0.132:
cs_c :=  0.2533/RS_FACTOR:
cs_d :=  0.349/RS_FACTOR:

thf := (ts, us, z) ->
  ((1.0 + z)/2.0)^(8.0/3.0)*(ts - us/8.0):

(* This is Equation (15) of Lee1988_785 *)
(* Note that gamma = (1 - z^2) *)
f_cs := (rs, z, xt, xs0, xs1, ts0, ts1, us0, us1) ->
  cs_a*(1.0 - z^2)/(1.0 + cs_d*rs) * (1.0 + 2.0*cs_b*exp(-cs_c*rs)*(
    thf(ts0, us0, z) + thf(ts1, us1, -z) - t_vw(z, xt, us0, us1)
  )):

f := (rs, z, xt, xs0, xs1, ts0, ts1, us0, us1) ->
  f_cs(rs, z, xt, xs0, xs1, ts0, ts1, us0, us1):
