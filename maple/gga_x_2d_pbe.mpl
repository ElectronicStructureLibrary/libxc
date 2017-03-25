(* type: work_gga_x *)

mkappa := 0.4604:
mmu    := 0.354546875:

f0 := s -> 1 + mkappa*(1 - mkappa/(mkappa + mmu*s^2)):
f  := x -> f0(X2S_2D*x):
