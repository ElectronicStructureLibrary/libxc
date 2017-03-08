(* type: work_mgga_c *)

$include "gga_c_scan_e0.mpl"
$include "mgga_x_scan.mpl"

scan_b1c := 0.0285764:
scan_b2c := 0.0889:
scan_b3c := 0.125541:
scan_eclda0 := rs -> -scan_b1c/(1.0 + scan_b2c*sqrt(rs) + scan_b3c*rs):

scan_chi_infty := 0.12802585262625815:
scan_g_infty := s -> 1.0/(1.0 + 4.0*scan_chi_infty*s^2)^(1.0/4.0):

(* in the paper it is 2.3631 *)
scan_G_cnst := 2.363:
scan_Gc := z -> (1.0 - scan_G_cnst*(2.0^(1.0/3.0) - 1.0)*f_zeta(z))*(1.0 - z^12):

scan_H0 := (rs, s) ->
  scan_b1c*log(1.0 + (exp(-scan_eclda0(rs)/scan_b1c) - 1.0)*(1.0 - scan_g_infty(s))):
scan_e0 := (rs, z, s) ->
  (scan_eclda0(rs) + scan_H0(rs, s))*scan_Gc(z):

scan_alpha := (z, xt, ts0, ts1) ->
  (t_total(z, ts0, ts1) - xt^2/8.0)/(K_FACTOR_C*t_total(z, 1.0, 1.0)):

(* set parameters of f_alpha *)
c1 := 0.64:
c2 := 1.5:
d  := 0.7:

f_scan := (rs, z, xt, xs0, xs1, ts0, ts1) ->
  f_pbe(rs, z, xt, xs0, xs1) + f_alpha(scan_alpha(z, xt, ts0, ts1))*(
    + scan_e0(rs, z, X2S*2.0^(1.0/3.0)*xt)
    - f_pbe(rs, z, xt, xs0, xs1)
  ):

f := (rs, z, xt, xs0, xs1, ts0, ts1, us0, us1) ->
  f_scan(rs, z, xt, xs0, xs1, ts0, ts1):
