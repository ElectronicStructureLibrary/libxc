(* type: work_mgga_x *)

a := 0.5389:
b := 3:

p := x -> X2S^2*x^2:
q := u -> X2S^2*u:

(* Equation (15) *)
fab := z -> piecewise( \
    z<=0, 0.0, \
    z>=a, 1.0, \
    (1.0 + exp(a/(a-z)))^b/(exp(a/z) + exp(a/(a-z)))^b \
):

(* Equation (7) *)
mDelta := (x, u) -> 8.0*q(u)^2/81.0 - p(x)*q(u)/9.0 + 8.0*p(x)^2/243.0:

f_W    := x -> 5.0*p(x)/3.0:

(* Equation (8) *)
f_GE4  := (x, u) -> 1.0 + 5.0*p(x)/27.0 + 20.0*q(u)/9.0 + mDelta(x, u):

(* Equation (11) *)
f_GE4_M := (x, u) -> f_GE4(x, u)/sqrt(1.0 + mDelta(x, u)^2/(1.0 + f_W(x))^2):

(* Equation (17) *)
f   := (rs, x, t, u) -> f_W(x) + (f_GE4_M(x, u) - f_W(x))*fab(f_GE4_M(x, u) - f_W(x)):
