(* type: work_gga_c *)
(* prefix:
  gga_c_lyp_params *params;

  assert(p->params != NULL);
  params = (gga_c_lyp_params * )(p->params);
*)

Cf := 3.0/10.0 * (3.0*Pi^2)^(2.0/3.0):

omega := rs -> params_a_B*exp(-params_a_c*rs)/(1.0 + params_a_d*rs):
delta := rs -> (params_a_c + params_a_d/(1.0 + params_a_d*rs))*rs:

aux6 := 1.0/2.0^(8.0/3.0):
aux4 := aux6/4.0:
aux5 := aux4/(9.0*2.0):

t1 := (rs, z) ->
  -(1.0 - z^2)/(1.0 + params_a_d*rs):
t2 := (rs, z, xt) ->
  -xt^2*((1.0 - z^2)*(47.0 - 7.0*delta(rs))/(4.0*18.0) - 2.0/3.0):
t3 := (z) ->
  -Cf/2.0*(1.0 - z^2)*((1.0 + z)^(8.0/3.0) + (1.0 - z)^(8.0/3.0)):
t4 := (rs, z, xs0, xs1) ->
  aux4*(1.0 - z^2)*(5.0/2.0 - delta(rs)/18.0)*(xs0^2*(1.0 + z)^(8.0/3.0) + xs1^2*(1.0 - z)^(8.0/3.0)):
t5 := (rs, z, xs0, xs1) ->
  aux5*(1.0 - z^2)*(delta(rs) - 11.0)*(xs0^2*(1.0 + z)^(11.0/3.0) + xs1^2*(1.0 - z)^(11.0/3.0)):
t6 := (z, xs0, xs1) ->
  -aux6*(2.0/3.0*(xs0^2*(1.0 + z)^(8.0/3.0) + xs1^2*(1.0 - z)^(8.0/3.0))
  -(1.0 + z)^2*xs1^2*(1.0 - z)^(8.0/3.0)/4.0 - (1.0 - z)^2*xs0^2*(1.0 + z)^(8.0/3.0)/4.0):

f_lyp := (rs, z, xt, xs0, xs1) -> params_a_A*(t1(rs/RS_FACTOR, z) + omega(rs/RS_FACTOR)*(
  + t2(rs/RS_FACTOR, z, xt) + t3(z) + t4(rs/RS_FACTOR, z, xs0, xs1)
  + t5(rs/RS_FACTOR, z, xs0, xs1) + t6(z, xs0, xs1)
)):

f  := (rs, z, xt, xs0, xs1) -> f_lyp(rs, z, xt, xs0, xs1):

