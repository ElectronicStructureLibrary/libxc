(* type: work_gga_c *)

$define lda_c_vbh_params
$include "lda_c_hl.mpl"

malpha := (4.0*Pi/3.0)^(1.0/6.0):

(* 4.28e-3/2.0, where the 2 comes from the covertion from Ryd. to Hartree *)
a1     := Pi/(16.0*(3.0*Pi^2)^(4.0/3.0)):
a2     := -0.262:
a3     := -7.0/(9.0*2.0^(5.0/3.0)):

t1 := (rs, z, xt) -> 
  + xt^2*exp(a2*xt/(malpha*sqrt(rs))) 
  * sqrt(2.0)/sqrt((1.0 + z)^(5.0/3.0) + (1.0 - z)^(5.0/3.0)):
t2 := (z, xs_0_, xs_1_) -> 
  + a3*(xs_0_^2*(1.0 + z)^(4.0/3.0) + xs_1_^2*(1.0 - z)^(4.0/3.0)):

f := (rs, z, xt, xs_0_, xs_1_) -> 
  + f_hl(rs, z)
  + a1/(malpha^2*rs)*(t1(rs, z, xt) + t2(z, xs_0_, xs_1_)):

