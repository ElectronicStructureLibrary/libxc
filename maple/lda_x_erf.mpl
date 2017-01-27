(* type: work_lda *)

params_a_alpha := 1.0:
$include "lda_x.mpl"

(* attenuation function *)
aux1_erf := a -> sqrt(Pi)*erf(1.0/(2.0*a)):
aux2_erf := a -> exp(-1.0/(4.0*a^2)) - 1.0:
aux3_erf := a -> 2.0*a^2*aux2_erf(a) + 0.5:
attenuation_erf := a -> 1.0 - 8.0/3.0*(aux1_erf(a) + 2.0*a*(aux2_erf(a) - aux3_erf(a))):

a_cnst := (4.0/(9.0*Pi))^(1.0/3.0)*p_a_cam_omega/2.0:

f    := (rs, z) -> ax*(
  + (1.0 + z)^(4.0/3.0)*attenuation_erf(a_cnst*rs/(1.0 + z)^(1/3))
  + (1.0 - z)^(4.0/3.0)*attenuation_erf(a_cnst*rs/(1.0 - z)^(1/3))
):