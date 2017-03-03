(* type: work_mgga_c *)

$define lda_c_pw_params
$include "lda_c_pw.mpl"

$define lda_x_params
$include "lda_x.mpl"

cc06_cnst  := (3.0/(4.0*Pi))^(2.0/3.0):

cc06_alpha := -0.0007:
cc06_beta  :=  0.0080*cc06_cnst:
cc06_gamma :=  0.026 *cc06_cnst:

f_cc06 := (rs, z, us0, us1) ->
  (f_lda_x(rs, z) + f_pw(rs, z))*(1.0 +
    (cc06_alpha + cc06_beta*u_total(z, us0, us1))/(1.0 + cc06_gamma*u_total(z, us0, us1))
  ):

f := (rs, z, xt, xs0, xs1, ts0, ts1, us0, us1) ->
  f_cc06(rs, z, us0, us1):
