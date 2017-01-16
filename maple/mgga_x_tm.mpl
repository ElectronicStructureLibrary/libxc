(* type: work_mgga_x *)

mlambda := 0.6866:
mbeta   := 79.873:

(* below Equation (6) *)
y  := p -> (2.0*mlambda - 1.0)^2 * p:

(* Equation (7) *)
f0 := x -> (1.0 + 10.0*70.0*y(X2S^2*x^2)/27.0 + mbeta*y(X2S^2*x^2)^2)^(1.0/10.0):

R  := (x, t) -> 1.0 + 595.0*(2.0*mlambda - 1.0)^2 * X2S^2*x^2/54.0 \
   - (t - 3.0*(mlambda^2 - mlambda + 0.5)*(t - K_FACTOR_C - x^2/72.0))/K_FACTOR_C:

fx_DME := (x, t) -> 1.0/f0(x)^2 + 7.0*R(x, t)/(9.0*f0(x)^4):

malpha := (x, t) -> (t - x^2/8.0)/K_FACTOR_C:
qtilde := (x, t) -> 9.0/20.0*(malpha(x, t) - 1.0) + 2.0*X2S^2*x^2/3.0:

fx_SC  := (x, t) -> (1.0 + 10.0*( \
       + (MU_GE + 50.0*X2S^2*x^2/729.0)*X2S^2*x^2 \
       + 146.0*qtilde(x, t)^2/2025.0 \
       - 73.0*qtilde(x,t)/405.0*(3.0*K_FACTOR_C/(5.0*t))*(1.0 - K_FACTOR_C/t)) \
       )^(1.0/10.0):

w := t-> (K_FACTOR_C^2/t^2 + 3.0*K_FACTOR_C^3/t^3)/(1.0 + K_FACTOR_C^3/t^3)^2:

f := (rs, x, t, u) -> w(t)*fx_DME(x, t) + (1.0 - w(t))*fx_SC(x, t):
