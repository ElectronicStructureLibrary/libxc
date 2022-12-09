(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)

herman_beta := 0.00314:

(*
  If I got the algebra right, in our notation the (unpolarized) LDA energy per atom is

    e_xc^LDA = (3/4) alpha V_XS

  where V_XS is defined in Eq. (7) of the IJQC, and alpha = 2/3.

  The energy expression is Eq. (8) of the IJQC article, with G(rho)
  defined in Eq. (9). The exponential damping term has been defined in
  Eq. (A11) of the Appendix, and the (1+beta*x^2) term in Eq. (A12).

*)

herman_f0 := (x, u) -> 1 + 2*herman_beta * (4/3*x^2 - 2*u)*(1+herman_beta*x^2)*exp(-herman_beta*x^2):

(* we have to take care of the spin factors coming out of the spin sum-rule *)

herman_f := (x, u, t) -> herman_f0(x/2^(1/3), u/2^(2/3)):

f := (rs, z, xt, xs0, xs1, u0, u1, t0, t1) ->
  mgga_exchange(herman_f, rs, z, xs0, xs1, u0, u1, t0, t1):
