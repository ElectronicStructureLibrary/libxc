(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: work_mgga_x *)

br89_a := [0, 1, 0, -2, 0, 1]:
fw := t -> mgga_series_w(br89_a, 6, t):

br89_gamma := 0.8:
br89_Q := (x, t, u) -> (u - 4*t*br89_gamma*Fermi_D(x, t))/6:
br89_y := (x, t, u) -> 2*Pi^(2/3)/(3*br89_Q(x, t, u)):

(* lower piece *)
pgk_a1 := 1.5255251812009530:
pgk_a2 := 0.4576575543602858:
pgk_a3 := 0.4292036732051034:

pgk_b  := [0.4771976183772063, -1.7799813494556270, 3.8433841862302150,
       -9.5912050880518490, 2.1730180285916720, -30.425133851603660]:

pgk_c  := [0.7566445420735584, -2.6363977871370960, 5.4745159964232880,
       -12.657308127108290, 4.1250584725121360, -30.425133957163840]:

pgk_d  := [0.00004435009886795587, 0.58128653604457910, 66.742764515940610,
       434.26780897229770, 824.7765766052239000, 1657.9652731582120]:

pgk_e  := [0.00003347285060926091, 0.47917931023971350, 62.392268338574240,
       463.14816427938120, 785.2360350104029000, 1657.962968223273000000]:

pgk_UB := 2.085749716493756:

pgk_x_lower := y -> (-arctan(pgk_a1*y + pgk_a2) + pgk_a3) *
            add(pgk_c[i]*y^(i-1), i=1..6)/add(pgk_b[i]*y^(i-1), i=1..6):

pgk_x_upper := y -> (arccsch(pgk_UB*y) + 2) *
            add(pgk_d[i]*y^(i-1), i=1..6)/add(pgk_e[i]*y^(i-1), i=1..6):
            
pgk_x := (x, t, u) -> convert(piecewise(
  br89_y(x, t, u) <= 0, pgk_x_lower(br89_y(x, t, u)),
  pgk_x_upper(br89_y(x, t, u))), 'Heaviside'):

br89_v_full := x -> exp(x/3.0)*(1.0 - exp(-x)*(1.0 + x/2.0))/x:
br89_v_expa := x -> 1/2 + x/6 - x^2/18:

br89_v := (x, t, u) ->
       -2*Pi^(1/3)/X_FACTOR_C * br89_v_full(pgk_x(x, t, u)):

f   := (rs, x, t, u) -> -br89_v(x, t, u)/2:
