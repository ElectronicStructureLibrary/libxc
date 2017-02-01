(* type: work_gga_c *)
(* prefix:
  gga_c_pbe_params *params;

  assert(p->params != NULL);
  params = (gga_c_pbe_params * )(p->params);
*)

$define lda_c_pw_params
$define lda_c_pw_modified_params
$include "lda_c_pw.mpl"

$ifdef gga_c_pbe_params
params_a_beta  := 0.06672455060314922:
params_a_gamma := (1.0 - log(2.0))/Pi^2:
params_a_BB    := 1.0:
$endif

mgamma := params_a_gamma:
mbeta  := (rs, t) -> params_a_beta:
BB     := params_a_BB:

mphi := z -> ((1.0 + z)^(2.0/3.0) + (1.0 - z)^(2.0/3.0))/2.0:
tt   := (rs, z, xt) -> xt/(4.0*2^(1.0/3.0)*mphi(z)*sqrt(rs)):
tp   := (rs, z, xt) -> tt(rs, z, xt):

(* Equation (8) *)
A := (rs, z, t) ->
  mbeta(rs, t)/(mgamma*(exp(-f_pw(rs, z)/(mgamma*mphi(z)^3)) - 1.0)):

(* Equation (7) *)
f1 := (rs, z, t) -> t^2 + BB*A(rs, z, t)*t^4:
f2 := (rs, z, t) -> mbeta(rs, t)*f1(rs, z, t)/(mgamma*(1.0 + A(rs, z, t)*f1(rs, z, t))):

fH := (rs, z, t) -> mgamma*mphi(z)^3*log(1.0 + f2(rs, z, t)):

f  := (rs, z, xt, xs_0_, xs_1_) ->
  f_pw(rs, z) + fH(rs, z, tp(rs, z, xt)):
