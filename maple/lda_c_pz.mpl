(* type: work_lda *)
(* prefix:
  lda_c_pz_params *params;
 
  assert(p->params != NULL);
  params = (lda_c_pz_params * )(p->params);
*)

(* Equation C3 *)
ec_low  := (i, rs) -> params_a_gamma[i] / \
        (1.0 + params_a_beta1[i]*sqrt(rs) + params_a_beta2[i]*rs):

(* Equation [1].C5 *)
ec_high := (i, rs) -> params_a_a[i]*log(rs) + params_a_b[i] \
        + params_a_c[i]*rs*log(rs) + params_a_d[i]*rs:

(* This is a little tricky as ec is discontinuous at rs=1, and therefore *)
(* it is not differentiable (and maple knows it). As a workaround, we    *)
(* write the function in terms of the Heaviside function, and handle the *)
(* Dirac functions that come out in the derivatives                      *)
ec := (i, x) -> convert(piecewise(x >= 1.0, ec_low(i, x), ec_high(i, x)), 'Heaviside'):

fzeta := z -> ((1.0 + z)^(4.0/3.0) + (1.0 - z)^(4.0/3.0) - 2.0)/(2.0^(4.0/3.0) - 2.0):

f := (rs, zeta) -> \
 ec(1, rs) + (ec(2, rs) - ec(1, rs))*fzeta(zeta):
