(* type: work_gga_x *)

mkappa := KAPPA_PBE:
mmu    := 0.00361218645365094697:
malpha := 0.52:

f  := x -> 1 + mkappa*(1 - (1 + mmu*x^2/(malpha*mkappa))^(-malpha)):
