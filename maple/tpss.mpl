tpss_csi2 := (z, xt, xs0, xs1) ->
  (1.0 - z^2)*(t_total(z, xs0^2, xs1^2) - xt^2)/(2.0*(3.0*Pi^2)^(1.0/3.0))^2:

tpss_C0 := (cc, z, xt, xs0, xs1) ->
  + add(cc[i]*z^(2*(i-1)), i=1..4)
  / (1.0 + tpss_csi2(z, xt, xs0, xs1)*((1.0 + z)^(-4.0/3.0) + (1.0 - z)^(-4.0/3.0))/2.0)^4:

tpss_aux := (z, xt, ts0, ts1) ->
  xt^2/(8.0*t_total(z, ts0, ts1)):

tpss_perp := (f_gga, rs, z, xt, xs0, xs1, ts0, ts1) ->
  (1.0 + tpss_C0(params_a_C0_c, z, xt, xs0, xs1)*tpss_aux(z, xt, ts0, ts1)^2)
  * f_gga(rs, z, xt, xs0, xs1):

tpss_par  := (f_gga, rs, z, xt, xs0, xs1, ts0, ts1) ->
  - (1.0 + tpss_C0(params_a_C0_c, z, xt, xs0, xs1))*tpss_aux(z, xt, ts0, ts1)^2*(
    + m_max(gga_stoll_par(f_gga, rs,  z, xs0,  1.0), f_gga(rs, z, xt, xs0, xs1))
    + m_max(gga_stoll_par(f_gga, rs, -z, xs1, -1.0), f_gga(rs, z, xt, xs0, xs1))
  ):

f_tpss0 := (f_gga, rs, z, xt, xs0, xs1, ts0, ts1) ->
  + tpss_par (f_gga, rs, z, xt, xs0, xs1, ts0, ts1)
  + tpss_perp(f_gga, rs, z, xt, xs0, xs1, ts0, ts1):

f_tpss := (f_gga, rs, z, xt, xs0, xs1, ts0, ts1) ->
  + f_tpss0(f_gga, rs, z, xt, xs0, xs1, ts0, ts1)
  * (1.0 + params_a_d*f_tpss0(f_gga, rs, z, xt, xs0, xs1, ts0, ts1)*tpss_aux(z, xt, ts0, ts1)^3)
:
