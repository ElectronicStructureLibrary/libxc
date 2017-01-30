(* type: work_gga_c *)

$define lda_c_pz_params
$include "lda_c_pz.mpl"

malpha := 0.023266:
mbeta  := 7.389e-6:
mgamma := 8.723:
mdelta := 0.472:
aa     := 0.001667:
bb     := 0.002568:
ftilde := 1.745*0.11:

(* Equation (4) *)
DD := z  -> sqrt((1.0 + z)^(5.0/3.0) + (1.0 - z)^(5.0/3.0))/sqrt(2.0):

(* Equation (6) *)
CC := rs -> 
  + aa 
  + (bb + malpha*rs + mbeta*rs^2)/(1.0 + mgamma*rs + mdelta*rs^2 + 1.0e4*mbeta*rs^3):
CCinf := aa + bb:

(* Equation (9) *)
x1   := (rs, xt) -> xt/sqrt(rs/RS_FACTOR):
mPhi := (rs, xt) -> ftilde*(CCinf/CC(rs))*x1(rs, xt):

(* Equation (8) *)
H := (rs, z, xt) -> x1(rs, xt)^2*exp(-mPhi(rs, xt))*CC(rs)/DD(z):

f := (rs, z, xt, xs_0_, xs_1_) ->
  f_pz(rs, z) + H(rs, z, xt):