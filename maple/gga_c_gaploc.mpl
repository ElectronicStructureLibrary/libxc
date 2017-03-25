(* type: work_gga_c *)

$include "gga_c_gapc.mpl"

(* The two parameters were fixed by fitting to the exact correlation
energy per particle of the He atom *)
gaploc_b  := 14.709046:
gaploc_a1 :=  6.54613 + 2:

gaploc_alpha := t  -> (gaploc_a1 + 3*t^3)/(1 + t^3):
gaploc_s     := xt -> X2S*xt*2^(1/3):

(* override definition of gap_C *)
(* The pre-factor is completely different from the one in
   Equation (7). I used the "gfac" from the original code. *)
gap_G := (rs, z, xt, par) -> (9*Pi/4)^(2/3)/2.0
  * gaploc_s(xt)^(gaploc_alpha(gap_t(rs, z, xt)))/rs^2
  * (gaploc_b + gaploc_s(xt)^2)/(1 + gaploc_s(xt)^(gaploc_alpha(gap_t(rs, z, xt)))):

