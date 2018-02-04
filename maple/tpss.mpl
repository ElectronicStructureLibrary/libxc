(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

tpss_csi2 := (z, xt, xs0, xs1) ->
  (1 - z^2)*(t_total(z, xs0^2, xs1^2) - xt^2)/(2*(3*Pi^2)^(1/3))^2:

tpss_C0 := (cc, z, xt, xs0, xs1) ->
  + add(cc[i]*z^(2*(i-1)), i=1..4)
  / (1 + tpss_csi2(z, xt, xs0, xs1)*((1 + z)^(-4/3) + (1 - z)^(-4/3))/2)^4:

tpss_aux := (z, xt, ts0, ts1) ->
  xt^2/(8*t_total(z, ts0, ts1)):

tpss_perp := (f_gga, rs, z, xt, xs0, xs1, ts0, ts1) ->
  (1 + tpss_C0(params_a_C0_c, z, xt, xs0, xs1)*tpss_aux(z, xt, ts0, ts1)^2)
  * f_gga(rs, z, xt, xs0, xs1):

tpss_par  := (f_gga, rs, z, xt, xs0, xs1, ts0, ts1) ->
  - (1 + tpss_C0(params_a_C0_c, z, xt, xs0, xs1))*tpss_aux(z, xt, ts0, ts1)^2*(
    + m_max(gga_stoll_par(f_gga, rs,  z, xs0,  1), f_gga(rs, z, xt, xs0, xs1))
    + m_max(gga_stoll_par(f_gga, rs, -z, xs1, -1), f_gga(rs, z, xt, xs0, xs1))
  ):

f_tpss0 := (f_gga, rs, z, xt, xs0, xs1, ts0, ts1) ->
  + tpss_par (f_gga, rs, z, xt, xs0, xs1, ts0, ts1)
  + tpss_perp(f_gga, rs, z, xt, xs0, xs1, ts0, ts1):

f_tpss := (f_gga, rs, z, xt, xs0, xs1, ts0, ts1) ->
  + f_tpss0(f_gga, rs, z, xt, xs0, xs1, ts0, ts1)
  * (1 + params_a_d*f_tpss0(f_gga, rs, z, xt, xs0, xs1, ts0, ts1)*tpss_aux(z, xt, ts0, ts1)^3)
:
