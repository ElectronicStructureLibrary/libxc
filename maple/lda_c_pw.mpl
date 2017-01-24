(* type: work_lda *)
(* prefix:
  lda_c_pw_params *params;

  assert(p->params != NULL);
  params = (lda_c_pw_params * )(p->params);
*)

(* Equation (10) *)
g_aux := (k, rs) -> params_a_beta1[k]*sqrt(rs) + params_a_beta2[k]*rs
  + params_a_beta3[k]*rs^1.5 + params_a_beta4[k]*rs^(params_a_pp[k] + 1.0):
g     := (k, rs) -> -2.0*params_a_a[k]*(1.0 + params_a_alpha1[k]*rs)
  * log(1.0 +  1.0/(2.0*params_a_a[k]*g_aux(k, rs))):

(* Equation (8) *)
(* Attention, the function g parametrizes -alpha *)
f_PW := (rs, zeta) ->
  g(1, rs) + zeta^4*f_zeta(zeta)*(g(2, rs) - g(1, rs) + g(3, rs)/params_a_fz20)
  - f_zeta(zeta)*g(3, rs)/params_a_fz20:

f := (rs, zeta) -> f_PW(rs, zeta):
