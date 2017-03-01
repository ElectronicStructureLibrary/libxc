
(* The B97 function g *)
b97_g := (gamma, cc, x) -> add(cc[i]*(gamma*x^2/(1.0 + gamma*x^2))^(i-1), i=1..5):

(* This is the stoll decomposition in our language *)
lda_stoll_par  := (lda_func, rs, z, spin) -> 
  lda_func(rs*(2.0/(1.0 + z))^(1.0/3.0), spin)*(1.0 + z)/2.0:

lda_stoll_perp := (lda_func, rs, z) ->
  + lda_func(rs, z) 
  - lda_stoll_par(lda_func, rs,  z,  1.0) 
  - lda_stoll_par(lda_func, rs, -z, -1.0):

(* The parallel and perpendicular components of the energy *)
b97_fpar  := (lda_func, mgamma, cc, rs, z, xs_0_, xs_1_) ->
  + lda_stoll_par(lda_func, rs,  z,  1.0) * b97_g(mgamma, cc, xs_0_) 
  + lda_stoll_par(lda_func, rs, -z, -1.0) * b97_g(mgamma, cc, xs_1_):

b97_fperp := (lda_func, mgamma, cc, rs, z, xs_0_, xs_1_) ->
  lda_stoll_perp(lda_func, rs, z) * b97_g(mgamma, cc, sqrt(xs_0_^2 + xs_1_^2)/sqrt(2.0)):

f_b97 := (lda_func, gamma_ss, cc_ss, gamma_ab, cc_ab, rs, z, xs_0_, xs_1_) ->
  + b97_fpar (lda_func, gamma_ss, cc_ss, rs, z, xs_0_, xs_1_)
  + b97_fperp(lda_func, gamma_ab, cc_ab, rs, z, xs_0_, xs_1_):