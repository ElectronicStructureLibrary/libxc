(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)

$define gga_x_b88_params
$include "gga_x_b88.mpl"

cab := 0.63:
css := 0.96:

b88_css := (rs, z, xs, ts) ->
  - 0.01 * (1 + z)/2 * n_spin(rs, z)^(5/3) * 2*ts * Fermi_D(xs, ts)
  * b88_zss(css, b88_f, rs, z, xs)^4 * (
    1 - 2*log(1 + b88_zss(css, b88_f, rs, z, xs)/2)
      / b88_zss(css, b88_f, rs, z, xs)
    ):

if evalb(Polarization = "ferr") then
  b88_par := (rs, z, xs0, xs1, ts0, ts1) ->
    + b88_css(rs,  1, xs0, ts0):

  b88_cab := (rs, z, xs0, xs1) -> 0:
else 
  b88_par := (rs, z, xs0, xs1, ts0, ts1) ->
    + b88_css(rs,  z, xs0, ts0)
    + b88_css(rs, -z, xs1, ts1):

  b88_cab := (rs, z, xs0, xs1) ->
    - 0.8 * (1 - z^2)/4 * n_total(rs)
    * b88_zab(cab, b88_f, rs, z, xs0, xs1)^2 * (
      1 - log(1 + b88_zab(cab, b88_f, rs, z, xs0, xs1))
        / b88_zab(cab, b88_f, rs, z, xs0, xs1)
      ):
end if:

b88_c_f := (rs, z, xs0, xs1, ts0, ts1) ->
  + b88_cab(rs,  z, xs0, xs1)
  + b88_par(rs,  z, xs0, xs1, ts0, ts1):

f := (rs, z, xt, xs0, xs1, us0, us1, ts0, ts1) ->
  b88_c_f(rs, z, xs0, xs1, ts0, ts1):
