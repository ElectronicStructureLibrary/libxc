(* Equation (7) from the paper *)

p     := x -> X2S^2*x^2:
z     := (x, t) -> x^2/(8.0*t):

alpha := (x, t) -> (t - x^2/8)/K_FACTOR_C:
qb    := (x, t) -> \
      9.0/20.0 * (alpha(x, t) - 1.0)/sqrt(1.0 + params_a_b*alpha(x, t)*(alpha(x, t) - 1.0)) \
      + 2.0*p(x)/3.0:

(* Equation (10) in all its glory *)
fxnum := (x, t) -> \
      + (MU_GE + params_a_c*z(x, t)^ff(z(x, t))/(1.0 + z(x, t)^2)^2)*p(x) \
      + 146.0/2025.0 * qb(x, t)^2 \
      - 73.0/405.0 * qb(x, t) * sqrt(0.5*(9.0/25.0*z(x, t)^2 + p(x)^2)) \
      + MU_GE^2/mkappa(x, t) * p(x)^2 \
      + 2.0*sqrt(params_a_e)*MU_GE*9.0/25.0*z(x, t)^2 \
      + params_a_e*params_a_mu*p(x)^3:

fxden := x -> \
      (1.0 + sqrt(params_a_e)*p(x))^2:

fx    := (x, t) -> fxnum(x, t)/fxden(x):
