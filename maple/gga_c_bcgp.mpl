(* type: work_gga_c *)

cac  := 1.467:
mtau := 4.5:
bcgp_pt := t -> t*sqrt((mtau + t)/(mtau + cac*t)):

$define gga_c_pbe_params
$include "gga_c_pbe.mpl"

(* override definition of tp *)
tp := (rs, z, xt) -> bcgp_pt(tt(rs, z, xt)):

(* override definition of A *)
A := (rs, z) -> 1:
