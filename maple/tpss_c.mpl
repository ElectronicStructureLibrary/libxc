(*
 Copyright (C) 2017 M.A.L. Marques
               2019 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* Equation 28 squared with Equation 25 built in *)
tpss_csi2 := (z, xt, xs0, xs1) ->
  (1 - z^2)*(t_total(z, xs0^2, xs1^2) - xt^2)/(2*(3*Pi^2)^(1/3))^2:

(* Equation 33 *)
tpss_C00 := (cc, z) ->
  + add(cc[i]*z^(2*(i-1)), i=1..4):

(* Equation 34 *)
tpss_C0 := (cc, z, xt, xs0, xs1) ->
  + tpss_C00(cc, z) / (1 + tpss_csi2(z, xt, xs0, xs1)*((1 + z)^(-4/3) + (1 - z)^(-4/3))/2)^4:

(* Equation 11, with tau_W from Equation 12 *)
tpss_aux := (z, xt, ts0, ts1) ->
  m_min(xt^2/(8*t_total(z, ts0, ts1)), 1):

if evalb(Polarization = "ferr") then
  tpss_par  := (f_gga, rs, z, xt, xs0, xs1, ts0, ts1) ->
    - tpss_aux(z, xt, ts0, 0)^2*f_gga(rs, 1, xt, xs0, 0):

  tpss_perp := (f_gga, rs, z, xt, xs0, xs1, ts0, ts1) -> 0:

else
  tpss_par  := (f_gga, rs, z, xt, xs0, xs1, ts0, ts1) ->
    - (1 + tpss_C0(params_a_C0_c, z, xt, xs0, xs1))*tpss_aux(z, xt, ts0, ts1)^2*(
      + m_max(f_gga(rs*(2/(1 + z))^(1/3),  1, xs0, xs0, 0), f_gga(rs, z, xt, xs0, xs1))*(1 + z)/2
      + m_max(f_gga(rs*(2/(1 - z))^(1/3), -1, xs1, 0, xs1), f_gga(rs, z, xt, xs0, xs1))*(1 - z)/2
    ):

  tpss_perp := (f_gga, rs, z, xt, xs0, xs1, ts0, ts1) ->
    (1 + tpss_C0(params_a_C0_c, z, xt, xs0, xs1)*tpss_aux(z, xt, ts0, ts1)^2)
    * f_gga(rs, z, xt, xs0, xs1):
end if:

tpss_f0 := (f_gga, rs, z, xt, xs0, xs1, ts0, ts1) ->
  + tpss_par (f_gga, rs, z, xt, xs0, xs1, ts0, ts1)
  + tpss_perp(f_gga, rs, z, xt, xs0, xs1, ts0, ts1):

(* Equation 24 *)
tpss_f := (f_gga, rs, z, xt, xs0, xs1, ts0, ts1) ->
  + tpss_f0(f_gga, rs, z, xt, xs0, xs1, ts0, ts1)
  * (1 + params_a_d*tpss_f0(f_gga, rs, z, xt, xs0, xs1, ts0, ts1)*tpss_aux(z, xt, ts0, ts1)^3):
