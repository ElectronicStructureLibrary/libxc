(* type: work_gga_c *)

$define gga_c_pbe_params
$include "gga_c_pbe.mpl"

(* in the paper we have beta_a = 0.066725 *)
beta_a := 0.066724550603149220:
beta_b := 0.1:
beta_c := 0.1778:

(* we redefine beta here *)
(* this is the Hu and Langreth expression *)
mbeta := (rs, t) -> beta_a*(1.0 + beta_b*rs)/(1.0 + beta_c*rs):