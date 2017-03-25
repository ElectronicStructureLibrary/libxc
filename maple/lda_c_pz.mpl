(* type: work_lda *)
(* prefix:
  lda_c_pz_params *params;

  assert(p->params != NULL);
  params = (lda_c_pz_params * )(p->params);
*)

$ifdef lda_c_pz_params
params_a_gamma := [-0.1423, -0.0843]:
params_a_beta1 := [ 1.0529,  1.3981]:
params_a_beta2 := [ 0.3334,  0.2611]:
params_a_a     := [ 0.0311,  0.01555]:
params_a_b     := [-0.048,  -0.0269]:
params_a_c     := [ 0.0020,  0.0007]:
params_a_d     := [-0.0116, -0.0048]:
$endif

(* Equation C3 *)
ec_low  := (i, rs) -> params_a_gamma[i] / \
        (1 + params_a_beta1[i]*sqrt(rs) + params_a_beta2[i]*rs):

(* Equation [1].C5 *)
ec_high := (i, rs) -> params_a_a[i]*log(rs) + params_a_b[i] \
        + params_a_c[i]*rs*log(rs) + params_a_d[i]*rs:

(* This is a little tricky as ec is discontinuous at rs=1, and therefore *)
(* it is not differentiable (and maple knows it). As a workaround, we    *)
(* write the function in terms of the Heaviside function, and handle the *)
(* Dirac functions that come out in the derivatives                      *)
ec := (i, x) -> convert(piecewise(x >= 1, ec_low(i, x), ec_high(i, x)), 'Heaviside'):

f_pz := (rs, zeta) -> \
 ec(1, rs) + (ec(2, rs) - ec(1, rs))*f_zeta(zeta):

f := (rs, zeta) -> f_pz(rs, zeta):