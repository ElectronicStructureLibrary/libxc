(* type: work_gga_c *)

$define lda_c_pw_params
$include "lda_c_pw.mpl"

pw91_C_c0  := 4.235e-3:
pw91_alpha := 0.09:
pw91_nu    := 16.0/Pi * (3.0*Pi^2)^(1.0/3.0):
pw91_beta  := pw91_nu*pw91_C_c0:

pw91_c1 := pw91_beta^2/(2.0*pw91_alpha):
pw91_c2 := 2.0*pw91_alpha/pw91_beta:

mphi := z -> ((1.0 + z)^(2.0/3.0) + (1.0 - z)^(2.0/3.0))/2.0:
tt   := (rs, z, xt) -> xt/(4.0*2.0^(1.0/3.0)*mphi(z)*sqrt(rs)):

(* Equation (14) *)
A := (rs, z) -> pw91_c2/(exp(-2.0*pw91_alpha*f_pw(rs, z)/(mphi(z)^3*pw91_beta^2)) - 1.0):

(* Equation (13) *)
H0 := (rs, z, t) -> pw91_c1*mphi(z)^3*log(1.0 
  + pw91_c2*(t^2 + A(rs, z)*t^4) / (1.0 + A(rs, z)*t^2 + A(rs, z)^2*t^4)
):

(* Pade parametrized form of C-xc found in
   M Rasolt & DJW Geldart, Phys. Rev. B 34, 1325 (1986)
*)
RS_a := [2.568, 23.266, 0.007389]:
RS_b := [1.0, 8.723, 0.472]:
RG_C_xc := rs -> (RS_a[1] + RS_a[2]*rs + RS_a[3]*rs^2)/(1000.0*(RS_b[1] + RS_b[2]*rs + RS_b[3]*rs^2)):

(* Equation (15) *)
C_xc0 := 2.568e-3:
C_x   := -0.001667:
h_a1  := -100.0 * 4.0/Pi * (4.0/(9.0*Pi))^(1.0/3.0):

H1 := (rs, z, t) -> pw91_nu * (RG_C_xc(rs) - C_xc0 - 3.0*C_x/7.0)
 * mphi(z)^3*t^2*exp(h_a1*rs*mphi(z)^4*t^2):

f  := (rs, z, xt, xs_0_, xs_1_) ->
  f_pw(rs, z) + H0(rs, z, tt(rs, z, xt)) + H1(rs, z, tt(rs, z, xt)):
