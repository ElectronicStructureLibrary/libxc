(* type: work_gga_c *)

C0 := (1.0 - log(2.0))/(2.0*Pi^2): (* Equation (9) *)
C1 := 4.0*C0/3.0:                  (* Equations (13), (28), (33) *)
C2 := RS_FACTOR:                   (* Equation (8) *)
C3 := C2/3.0:

big    := 1e4:
cutoff := 1.0e7:

(* Equation (39) *)
kssp0_k0 := 1.291551074:
kssp0_k1 := 0.349064173:
kssp0_r1 := 0.08327588:

kssp0 := rs -> convert(piecewise(rs > big, kssp0_k0 - kssp0_k1, 
  kssp0_k0 - kssp0_k1*(1.0 - exp(-kssp0_r1*rs^(4.0/5.0)))
), 'Heaviside'):

(* Equation (45) *)
fssp_A1 := 1.622118767:
fssp_A2 := 0.489958076:
fssp_A3 := 1.379021941:

fssp := (rs, gr) -> convert(piecewise(fssp_A2^2*gr^2 > big, 0.0,
  (1.0 + fssp_A1*gr + fssp_A2^2*gr^2)*exp(-fssp_A2^2*gr^2)/sqrt(1.0 + fssp_A3*gr/rs)
), 'Heaviside'):

(* Equation (34) *)
fa_a1 := 0.939016:
fa_a2 := 1.733170:

f_factor := rs -> exp(-rs^2/(fa_a1*sqrt(rs) + fa_a2*rs)^2):

(* Equation (40) *)
kss0_k0 :=  1.200801774:
kss0_k1 :=  0.859614445:
kss0_k2 := -0.812904345:
kss0_r1 :=  1.089338848:
kss0_r2 :=  0.655638823:

kss0 := (rs, gr) -> convert(piecewise(rs > big, kss0_k0 + kss0_k1 + kss0_k2,
 kss0_k0 + kss0_k1*(1.0 - exp(-kss0_r1*sqrt(rs))) + kss0_k2*(1.0 - exp(-kss0_r2*rs^(2.0/5.0)))
), 'Heaviside'):

fss_A4 := 4.946281353:
fss_A5 := 3.600612059:

fss := (rs, gr) -> convert(piecewise(fss_A4^2*gr^2 > big, 0.0,
  (1.0 + fss_A4^2*gr^2)*exp(-fss_A4^2*gr^2)/sqrt(1.0 + fss_A5*gr/rs)
), 'Heaviside'):

(* Equation (15) *)
eq15 := mu -> (3.0 + 2.0*(sqrt(mu) + mu))/(3.0 + 6.0*(sqrt(mu) + mu)):

f_eab := mu -> C0*(exp(mu)*Ei(-mu) + 2.0*eq15(mu)*(mu*exp(mu)*Ei(-mu) + 1.0)):

(* Equation (13) *)
mu_ba := (rsa, ga2) -> C1*rsa/(kssp0(rsa)^2*fssp(rsa, ga2)^2):
term1 := (rsa, ga2, z) -> piecewise(mu_ba(rsa, ga2) > cutoff, 0.0,
  f_eab(mu_ba(rsa, ga2))*(1.0 + z)/2.0
):

mu_aa := (rsa, ga2) -> C1*rsa/(kss0(rsa)^2*fss(rsa, ga2)^2):
term2 := (rsa, z, ga2) -> piecewise(mu_aa(rsa, ga2) > cutoff, 0.0,
  f_eab(mu_aa(rsa, ga2))*f_factor(rsa)*(1.0 + z)/2.0
):

f_ft97 := (rs, z, xs) ->
  + term1(rs*2.0^(1.0/3.0)/(1.0 + z)^(1.0/3.0), z, C3^2*xs^2)
  + term2(rs*2.0^(1.0/3.0)/(1.0 + z)^(1.0/3.0), z, C3^2*xs^2):

fa  := (rs, z, xt, xs_0_, xs_1_) -> 
  f_ft97(rs, z, xs_0_) + f_ft97(rs, -z, xs_1_):

f  := (rs, z, xt, xs_0_, xs_1_) -> 
  term2(rs*2.0^(1.0/3.0)/(1.0 + z)^(1.0/3.0), z, C3^2*xs_0_^2):

