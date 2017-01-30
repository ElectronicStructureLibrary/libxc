(* type: work_gga_c *)

$define lda_c_pw_modified_params
$include "lda_c_pw.mpl"

malpha := 2.804:
mgamma := 0.8098:

XX := s -> 1.0/(1.0 + malpha*s^2):
ff := s -> XX(s) + mgamma*(1.0 - XX(s)):

f := (rs, z, xt, xs_0_, xs_1_) -> f_pw(rs, z)*(
  + (1.0 + z)/2.0 * ff(X2S*xs_0_)
  + (1.0 - z)/2.0 * ff(X2S*xs_1_)
):