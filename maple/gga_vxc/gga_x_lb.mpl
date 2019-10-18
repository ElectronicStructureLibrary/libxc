(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_vxc *)
(* prefix:
  gga_x_lb_params *params;

  assert(p->params != NULL);
  params = (gga_x_lb_params * )(p->params);
*)

lb_f := (rs, z, x) -> params_a_alpha*lda_x_spin(rs, z) -
  my_piecewise3(x < 300,
              params_a_beta*x^2/(1 + 3*params_a_beta*x*arcsinh(params_a_gamma*x)),
              x/(3.0*log(2*params_a_gamma*x))):

f := (rs, z, xt, xs0, xs1) -> lb_f(rs, z, xs0):
