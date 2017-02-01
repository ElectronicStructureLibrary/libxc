(* type: work_gga_c *)

$define gga_c_pbe_params
$include "gga_c_pbe.mpl"

pbeloc_b0 := 0.0375:
pbeloc_a  := 0.08:

(* we redefine beta here *)
mbeta := (rs, t) -> pbeloc_b0 + pbeloc_a*t^2*(1.0 - exp(-rs^2)):