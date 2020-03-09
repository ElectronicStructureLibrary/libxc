(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: lda_exc *)
(* replace: "int1\(" -> "xc_integrate(func1, NULL, 0.0, " *)
(* replace: "int2\(" -> "xc_integrate(func2, NULL, 0.0, " *)
(* prefix:
  lda_x_1d_exponential_params *params;

  assert(p->params != NULL);
  params = (lda_x_1d_exponential_params * )(p->params);
*)

$define xc_dimensions_1d

`diff/int1` := proc(g, x) diff(g, x) * f_inter(g)   end proc:
`diff/int2` := proc(g, x) diff(g, x) * f_inter(g)*g end proc:

f_inter := x -> xc_E1_scaled(x^2):

x1d_R := (rs, z) -> Pi*params_a_beta/(2*rs):
x1d_f_spin := (rs, z) ->
  -((1 + z)*int1((1 + z)*x1d_R(rs, z)) - int2((1 + z)*x1d_R(rs, z))/x1d_R(rs, z))
    / (4.0*Pi*params_a_beta):

if evalb(Polarization = "ferr") then
  f := (rs, z) -> x1d_f_spin(rs, 1):
else
  f := (rs, z) -> x1d_f_spin(rs, z) + x1d_f_spin(rs, -z):
end if:
