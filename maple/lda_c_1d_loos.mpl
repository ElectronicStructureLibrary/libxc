(* type: work_lda *)

eta0 := -log(sqrt(2.0*Pi)) - 3.0/4.0:
eta1 := 0.359933:
eps0 := -Pi^2/360.0:
eps1 := 0.00714:

kappa := 0.3083:

c0 := kappa*eta0:
c1 := 4.0*kappa*eta0 + kappa*sqrt(kappa)*eta1:
c2 := 5.0*eps0 + eps1/kappa:
c3 := eps1:

t := rs -> (sqrt(1.0 + 4.0*kappa*rs) - 1.0)/(2.0*kappa*rs):

f := (rs, z) ->
  t(rs)^2*(c0*(1.0 - t(rs))^3 + c1*t(rs)*(1.0 - t(rs))^2 + c2*t(rs)^2*(1.0 - t(rs)) + c3*t(rs)^3):
