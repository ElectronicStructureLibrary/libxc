(* Equation (7) from the paper *)

alpha := (p, z) -> 5.0/3.0*p*(1.0 - z)/z:
qb    := (p, z) -> \
      9.0/20.0 * (alpha(p, z) - 1.0)/sqrt(1.0 + params_a_b*alpha(p, z)*(alpha(p, z) - 1.0)) \
      + 2.0*p/3.0:

(* Equation (10) in all its glory *)
fxnum := (p, z) -> \
      + (MU_GE + params_a_c*z^ff(z)/(1.0 + z^2)^2)*p \
      + 146.0/2025.0 * qb(p, z)^2 \
      - 73.0/405.0 * qb(p, z) * sqrt(0.5*(9.0/25.0*z^2 + p^2)) \
      + MU_GE^2/params_a_kappa * p^2 \
      + 2.0*sqrt(params_a_e)*MU_GE*9.0/25.0*z^2 \
      + params_a_e*params_a_mu*p^3:

fxden := p -> \
      (1.0 + sqrt(params_a_e)*p)^2:

fx    := (x, t) -> fxnum(X2S^2*x^2, x^2/(8.0*t))/fxden(X2S^2*x^2):
