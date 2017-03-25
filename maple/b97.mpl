
(* The B97 function g *)
b97_g := (gamma, cc, x) -> add(cc[i]*(gamma*x^2/(1 + gamma*x^2))^(i-1), i=1..5):

(* The parallel and perpendicular components of the energy *)
b97_fpar  := (lda_func, mgamma, cc, rs, z, xs0, xs1) ->
  + lda_stoll_par(lda_func, rs,  z,  1) * b97_g(mgamma, cc, xs0) 
  + lda_stoll_par(lda_func, rs, -z, -1) * b97_g(mgamma, cc, xs1):

b97_fperp := (lda_func, mgamma, cc, rs, z, xs0, xs1) ->
  lda_stoll_perp(lda_func, rs, z) * b97_g(mgamma, cc, sqrt(xs0^2 + xs1^2)/sqrt(2)):

f_b97 := (lda_func, gamma_ss, cc_ss, gamma_ab, cc_ab, rs, z, xs0, xs1) ->
  + b97_fpar (lda_func, gamma_ss, cc_ss, rs, z, xs0, xs1)
  + b97_fperp(lda_func, gamma_ab, cc_ab, rs, z, xs0, xs1):
