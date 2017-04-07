(* type: work_mgga_x *)

$include "mgga_x_gx.mpl"

mmu := 0.001015549:
pbe_gx := x -> 1/(1 + mmu*x^2):

f := (rs, x, t, u) ->
  f_gx_a(malpha(x, t)) * pbe_gx(x):
