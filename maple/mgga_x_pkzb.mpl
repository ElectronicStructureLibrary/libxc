(* type: work_mgga_x *)

a2 := 146.0/2025.0:
a3 := -73.0/405.0:
a4 := 0.131957187845257783631757384393: (* DD + 100.0/(81.0*81.0*kappa) *)

kappa := 0.804:

xx := (p, qt) -> MU_GE*p + a2*qt*qt + a3*qt*p + a4*p*p:

f := (rs, x, t, u) -> 1.0 + kappa - \
  kappa^2/(kappa + xx(X2S^2*x^2, 6.0*X2S^2*t - 9.0/20.0 - X2S^2*x^2/12.0)):
