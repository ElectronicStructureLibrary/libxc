(*
 Copyright (C) 2017 M.A.L. Marques
               2022 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)
(* prefix:
  hyb_mgga_xc_gas22_params *params;

  assert(p->params != NULL);
  params = (hyb_mgga_xc_gas22_params * )(p->params);
*)

(* the functional form is nearly identical to wb97mv *)
$include "hyb_mgga_xc_wb97mv.mpl"

(* these definitions must come after inserting the file *)
gas22_par_nx := 3:
(* More accurate value from the ipynb notebook *)
gas22_gamma_x := 0.003840616724010807:
gas22_par_x := [
    [  params_a_c_x[1], 0, 0],
    [  params_a_c_x[2], 0, 1],
    [  params_a_c_x[3], 1, 0]
]:

gas22_par_nss := 5:
(* More accurate value from the ipynb notebook *)
gas22_gamma_ss := 0.46914023462026644:
gas22_par_ss := [
    [  params_a_c_ss[1], 0, 1],
    [  params_a_c_ss[2], 1, 0],
    [  params_a_c_ss[3], 2, 0],
    [  params_a_c_ss[4], 0, 6],
    [  params_a_c_ss[5], 4, 6]
]:

gas22_par_nos := 5:
gas22_par_os := [
    [ params_a_c_os[1], 0, 0],
    [ params_a_c_os[2], 0, 2],
    [ params_a_c_os[3], 0, 6],
    [ params_a_c_os[4], 2/3, 6],
    [ params_a_c_os[5], 2/3, 2]
]:

(* the peculiarity of GAS22 is that it uses a slightly different expansion for
the anti-parallel part of the correlation in terms of xt and not of ux *)
b97mv_ux_os := (mgamma, x) -> x:
