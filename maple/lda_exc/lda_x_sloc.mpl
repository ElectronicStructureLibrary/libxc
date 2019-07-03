(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: lda_exc *)

# https://onlinelibrary.wiley.com/doi/full/10.1002/qua.25312

sloc_a := 1.67:
sloc_b := 0.3:

f_sloc := (rs, z) -> 
  -sloc_a/(2*(sloc_b + 1)) * n_total(rs)^sloc_b *
  ((1 + z)^(sloc_b + 1) + (1 - z)^(sloc_b + 1)):

f := (rs, z) -> f_sloc(rs, z):
