(*
 Copyright (C) 2018 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: work_gga_x *)

(* equation 1 *)
f_chachiyo := x -> (3*x^2 + Pi^2*log(x+1)) / ((3*x + Pi^2)*log(x+1)):

(* x has a different definition in the Chachiyo paper  *)
x_cha := x -> 2/9 * (Pi/3)**(1/3) * (2**(-1/3) * x):

f := x -> f_chachiyo(x_cha(x)):
