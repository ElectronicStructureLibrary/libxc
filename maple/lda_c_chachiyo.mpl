(* type: work_lda *)

ap := -0.01554535:
bp := 20.4562557:
af := ap/2.0:
bf := 27.4203609:

e0 := rs -> ap*log(1.0 + bp/rs + bp/rs^2):
e1 := rs -> af*log(1.0 + bf/rs + bf/rs^2):

f := (rs, zeta) -> e0(rs) + (e1(rs) - e0(rs))*f_zeta(zeta):
