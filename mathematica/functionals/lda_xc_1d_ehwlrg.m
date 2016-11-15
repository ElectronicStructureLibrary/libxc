(* type: work_lda *)
(* prefix:
  lda_xc_1d_ehwlrg_params *params;
 
  assert(p->params != NULL);
  params = (lda_xc_1d_ehwlrg_params * )(p->params);
*)

n[rs_] = 1/(2 rs)

f[rs_, zeta_] = (paramsa1 + paramsa2 n[rs] + paramsa3 n[rs]^2) n[rs]^paramsalpha
