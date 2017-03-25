(* type: work_gga_c *)
(* prefix:
  gga_c_sogga11_params *params;

  assert(p->params != NULL);
  params = (gga_c_sogga11_params * )(p->params);
*)

$define lda_c_pw_params
$define lda_c_pw_modified_params
$include "lda_c_pw.mpl"

mbeta  := 15.75592*0.004235: (* the usual value of 0.066726 *)
malpha := mbeta/(16*2^(2/3)):

yy := (rs, z, xt) -> -malpha*mphi(z)*xt^2/(rs*f_pw(rs, z)):

f0 := (rs, z, xt) -> 1 - 1/(1 + yy(rs, z, xt)):
f1 := (rs, z, xt) -> 1 - exp(-yy(rs, z, xt)):

t0 := (rs, z, xt) -> add(params_a_sogga11_a[i]*f0(rs, z, xt)^(i-1), i=1..6):
t1 := (rs, z, xt) -> add(params_a_sogga11_b[i]*f1(rs, z, xt)^(i-1), i=1..6):

f_sogga11 := (rs, z, xt, xs0, xs1) ->
  f_pw(rs, z)*(t0(rs, z, xt) + t1(rs, z, xt)):

f  := (rs, z, xt, xs0, xs1) -> f_sogga11(rs, z, xt, xs0, xs1):