(* type: mgga_exc *)
(* prefix:
*)

$define lda_c_pw_params
$include "lda_c_pw.mpl"

# the localized mass
ldms_ms := (z, ts0, ts1) -> K_FACTOR_C/(2**(2/3) * t_total(z, ts0, ts1)):

# the LDA-m* functional
ldms_f := (rs, z, ts0, ts1) ->
  ldms_ms(z, ts0, ts1) * f_pw(ldms_ms(z, ts0, ts1)*rs, z):

f := (rs, z, xt, xs0, xs1, us0, us1, ts0, ts1) ->
  ldms_f(rs, z, ts0, ts1):
