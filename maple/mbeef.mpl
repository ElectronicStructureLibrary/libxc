k  := 6.5124: (* PBEsol transformation *)
xi := p -> 2*p/(k + p) - 1:
xj := a -> - (1 - a^2)^3/(1 + a^3*(1 + a^3)):

with(orthopoly):

mgga_mbeef := (x, t) -> add(add(
  + coefs_mbeef[i][j]
  * P(j-1, xi(X2S^2*x^2))
  * P(i-1, xj((t - x^2/8)/K_FACTOR_C)),
i=1..n_mbeef), j=1..n_mbeef):
