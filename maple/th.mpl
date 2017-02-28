(* type: work_gga_c *)

XX := (z, xs) -> xs*((1.0 + z)/2.0)^(4.0/3.0):
YY := (z, xt, xs_0_, xs_1_) -> 2.0*(XX(z, xs_0_)^2 + XX(-z, xs_1_)^2) - xt^2:

f_th := (rs, z, xt, xs_0_, xs_1_) -> add(params_a_omega[i]
  * (n_spin(rs, z)^params_a_a[i] + n_spin(rs, -z)^params_a_a[i])
  * z^(2*params_a_b[i])
  * 0.5*(XX(z, xs_0_)^params_a_c[i] + XX(-z, xs_1_)^params_a_c[i])
  * YY(z, xt, xs_0_, xs_1_)^params_a_d[i], i=1..params_a_n)/n_total(rs):