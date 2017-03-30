(* type: work_gga_x *)

(* constants from text in section 2 *)
kappa_pbe := 0.814:
kappa_revpbe := 1.227:
mu := 0.219:

(* parameters from section 4 *)
alpha := 2*(3*Pi^2)^(1/3):
beta := 3:

(* eq 3 *)
fermi := s -> 1/(1+exp(-alpha*(s-beta))):
(* eq 5 *)
kappa := s -> (1-fermi(s))*kappa_revpbe + fermi(s)*kappa_pbe:
(* eq 4 *)
f0_pbetrans := s -> 1 + kappa(s)*(1 - kappa(s)/(kappa(s) + mu*s^2)):

f_pbetrans  := x -> f0_pbetrans(X2S*x):
f  := x -> f_pbetrans(x):

