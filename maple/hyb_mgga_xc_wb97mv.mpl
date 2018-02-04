(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: work_mgga_c *)

par_n := 6:

gamma_x := 0.004:
par_x := [
    [  0.85,  0, 0],
    [  1.007, 0, 1],
    [  0.259, 1, 0],
    [  0,   0, 0],
    [  0,   0, 0],
    [  0,   0, 0]
]:

gamma_ss := 0.2:
par_ss := [
    [  0.443,  0, 0],
    [ -1.437,  0, 4],
    [ -4.535,  1, 0],
    [ -3.39,   2, 0],
    [  4.278,  4, 3],
    [  0,    0, 0]
]:

gamma_os := 0.006:
par_os := [
    [  1.000,  0, 0],
    [  1.358,  1, 0],
    [  2.924,  2, 0],
    [ -8.812,  2, 1],
    [ -1.39,   6, 0],
    [  9.142,  6, 1]
]:

$include "lda_x_erf.mpl"
$include "b97mv.mpl"

f :=  (rs, z, xt, xs0, xs1, ts0, ts1, us0, us1) ->
  f_b97mv(f_lda_x_erf, rs, z, xt, xs0, xs1, ts0, ts1):
