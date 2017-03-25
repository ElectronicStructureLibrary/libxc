(* type: work_gga_c *)

$define lda_c_pw_params
$define lda_c_pw_modified_params
$include "lda_c_pw.mpl"

malpha := 2.804:
mgamma := 0.8098:

XX := s -> 1/(1 + malpha*s^2):
ff := s -> XX(s) + mgamma*(1 - XX(s)):

f := (rs, z, xt, xs0, xs1) -> f_pw(rs, z)*(
  + (1 + z)/2 * ff(X2S*xs0)
  + (1 - z)/2 * ff(X2S*xs1)
):