(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

mbeef_k  := 6.5124: (* PBEsol transformation *)
mbeef_xi := p -> 2*p/(mbeef_k + p) - 1:

(* expansion is better than 1e-13 for large a *)
mbeef_xj_threshold := 1e4:
mbeef_xj_normal := a -> - (1 - a^2)^3/(1 + a^3*(1 + a^3)):
mbeef_xj_large  := a -> 1 - 3/a^2 - 1/a^3:
mbeef_xj := a -> my_piecewise3(a > mbeef_xj_threshold,
  mbeef_xj_large (m_max(a, mbeef_xj_threshold)),
  mbeef_xj_normal(m_min(a, mbeef_xj_threshold))
):

with(orthopoly):

mbeef_expansion := (x, t) -> add(add(
  + mbeef_coefs[i][j]
  * P(j-1, mbeef_xi(X2S^2*x^2))
  * P(i-1, mbeef_xj((t - x^2/8)/K_FACTOR_C)),
i=1..mbeef_n), j=1..mbeef_n):
