(* type: work_mgga_x *)

mlambda := 0.6866:
mbeta   := 79.873:

(* below Equation (6) *)
y  := p -> (2*mlambda - 1)^2 * p:

(* Equation (7) *)
f0 := x -> (1 + 10*70*y(X2S^2*x^2)/27 + mbeta*y(X2S^2*x^2)^2)^(1/10):

R  := (x, t) -> 1 + 595*(2*mlambda - 1)^2 * X2S^2*x^2/54 \
   - (t - 3*(mlambda^2 - mlambda + 1/2)*(t - K_FACTOR_C - x^2/72))/K_FACTOR_C:

fx_DME := (x, t) -> 1/f0(x)^2 + 7*R(x, t)/(9*f0(x)^4):

malpha := (x, t) -> (t - x^2/8)/K_FACTOR_C:
qtilde := (x, t) -> 9/20*(malpha(x, t) - 1) + 2*X2S^2*x^2/3:

fx_SC  := (x, t) -> (1 + 10*( \
       + (MU_GE + 50*X2S^2*x^2/729)*X2S^2*x^2 \
       + 146*qtilde(x, t)^2/2025 \
       - 73*qtilde(x,t)/405*(3*K_FACTOR_C/(5*t))*(1 - K_FACTOR_C/t)) \
       )^(1/10):

w := t-> (K_FACTOR_C^2/t^2 + 3*K_FACTOR_C^3/t^3)/(1 + K_FACTOR_C^3/t^3)^2:

f := (rs, x, t, u) -> w(t)*fx_DME(x, t) + (1 - w(t))*fx_SC(x, t):
