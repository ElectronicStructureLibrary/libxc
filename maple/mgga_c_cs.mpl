(* type: work_mgga_c *)

cs_a := -0.04918:
cs_b :=  0.132:
cs_c :=  0.2533/RS_FACTOR:
cs_d :=  0.349/RS_FACTOR:

thf := (ts, us, z) ->
  ((1.0 + z)/2.0)^(8.0/3.0)*(ts - us/8.0):
tw  := (xt, z, us_0_, us_1_) ->
  xt^2/8.0 - ((1.0 + z)/2.0)^(5.0/3.0)*us_0_/8.0 - ((1.0 - z)/2.0)^(5.0/3.0)*us_1_/8.0:

(* This is Equation (15) of Lee1988_785 *)
(* Note that gamma = (1 - z^2) *)
f_cs := (rs, z, xt, xs_0_, xs_1_, ts_0_, ts_1_, us_0_, us_1_) ->
  cs_a*(1.0 - z^2)/(1.0 + cs_d*rs) * (1.0 + 2.0*cs_b*exp(-cs_c*rs)*(
    thf(ts_0_, us_0_, z) + thf(ts_1_, us_1_, -z) - tw(xt, z, us_0_, us_1_)
  )):

f := (rs, z, xt, xs_0_, xs_1_, ts_0_, ts_1_, us_0_, us_1_) ->
  f_cs(rs, z, xt, xs_0_, xs_1_, ts_0_, ts_1_, us_0_, us_1_):
