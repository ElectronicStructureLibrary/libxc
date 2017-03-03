(* type: work_gga_c *)
(* prefix:
  gga_x_n12_params *params;

  assert(p->params != NULL);
  params = (gga_x_n12_params * )(p->params);
*)

omega_x := 2.5:
gamma_x := 0.004:

rss := (rs, z) -> rs * (2.0/(1.0 + z))^(1.0/3.0):

vx := rs -> 1.0/(1.0 + (1.0/(RS_FACTOR*omega_x))*rs):
ux := x -> gamma_x*x^2/(1.0 + gamma_x*x^2):

FN12 := (rs, x) -> 
  + add(1.0*params_a_CC_0_[i+1]*ux(x)^i, i=0..3)
  + add(1.0*params_a_CC_1_[i+1]*ux(x)^i, i=0..3) * vx(rs)
  + add(1.0*params_a_CC_2_[i+1]*ux(x)^i, i=0..3) * vx(rs)^2
  + add(1.0*params_a_CC_3_[i+1]*ux(x)^i, i=0..3) * vx(rs)^3:

f_n12 := (rs, z, xt, xs0, xs1) -> -X_FACTOR_C*RS_FACTOR*(
  + (1.0 + z)*FN12(rss(rs,  z), xs0)/(2.0*rss(rs,  z))
  + (1.0 - z)*FN12(rss(rs, -z), xs1)/(2.0*rss(rs, -z))
):

f  := (rs, z, xt, xs0, xs1) -> f_n12(rs, z, xt, xs0, xs1):