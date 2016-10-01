(* type: work_lda *)
(* prefix:
  lda_x_params *params;
 
  assert(p->params != NULL);
  params = (lda_x_params * )(p->params);
*)

ax = -paramsalpha 3/4 (3/(2*Pi))^(2/3)

f[rs_, zeta_] = ax ((1 + zeta)^(4.0/3.0) + (1 - zeta)^(4.0/3.0)) / rs
