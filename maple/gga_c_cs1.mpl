(* type: work_gga_c *)

cs1_gamma :=  0.006: (* as in B88 *)
cs1_d     :=  0.349: (* as in CS  *)
cs1_C1    :=  0.018897:
cs1_C2    := -0.155240:
cs1_C3    :=  0.159068:
cs1_C4    := -0.007953:

(* Equation (24) *)
cs1_ess := (rs, z, xs) ->
  + (1.0 + z)/2.0 * 1.0/(1.0 + cs1_d*n_spin(rs, z)^(-1/3))
  * (cs1_C1 + cs1_C2*xs^4/(1.0 + cs1_gamma*xs^2)^2):

(* Equation (25) *)
cs1_eab := (rs, z, xt) ->
  - (1.0 - z^2)/4.0 * 1.0/(1.0 + cs1_d*n_total(rs)^(-1/3))
  * (cs1_C3 + cs1_C4*xt^4/(1.0 + cs1_gamma*xt^2)^2):

f_cs1 := (rs, z, xt, xs0, xs1) ->
  + cs1_eab(rs,  z, xt)
  + cs1_ess(rs,  z, xs0)
  + cs1_ess(rs, -z, xs1):

f := (rs, z, xt, xs0, xs1) ->
  f_cs1(rs, z, xt, xs0, xs1):
