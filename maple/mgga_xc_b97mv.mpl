(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: work_mgga_c *)

par_n := 5:

gamma_x := 0.004:
par_x := [
    [  1.000, 0, 0],
    [  1.308, 0, 1],
    [  1.901, 0, 2],
    [  0.416, 1, 0],
    [  3.070, 1, 1]
]:

gamma_ss := 0.2:
par_ss := [
    [  1.000, 0, 0],
    [ -1.855, 0, 2],
    [ -5.668, 1, 0],
    [-20.497, 3, 2],
    [-20.364, 4, 2]
]:

gamma_os := 0.006:
par_os := [
    [  1.000, 0, 0],
    [  1.573, 0, 1],
    [ -6.298, 0, 3],
    [  2.535, 1, 0],
    [ -6.427, 3, 2]
]:

$define lda_x_params
$include "lda_x.mpl"
$include "b97mv.mpl"

f :=  (rs, z, xt, xs0, xs1, ts0, ts1, us0, us1) ->
  f_b97mv(f_lda_x, rs, z, xt, xs0, xs1, ts0, ts1):
