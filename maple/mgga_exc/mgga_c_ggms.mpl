(* type: mgga_exc *)
(* prefix:
*)

$include "mgga_c_ldms.mpl"
$define gga_c_pbe_params
$include "gga_c_pbe.mpl"

# the GGA-m* functional
ggms_f := (rs, z, xt, ts0, ts1) -> ldms_ms(z, ts0, ts1) * (
  f_pw(ldms_ms(z, ts0, ts1)*rs, z) + fH(ldms_ms(z, ts0, ts1)*rs, z, tp(ldms_ms(z, ts0, ts1)*rs, z, xt))
):

f := (rs, z, xt, xs0, xs1, us0, us1, ts0, ts1) ->
  ggms_f(rs, z, xt, ts0, ts1):
