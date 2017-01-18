(* type: work_gga_x *)

mkappa := KAPPA_PBE:
mmu    := 0.00361218645365094697:
malpha := 0.52:

f  := x -> 1.0 + mkappa*(1.0 - (1.0 + mmu*x^2/(malpha*mkappa))^(-malpha)):
