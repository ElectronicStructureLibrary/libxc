(* type: work_gga_c *)
(* prefix:
  gga_x_hjs_params *params;
 
  assert(p->params != NULL);
  params = (gga_x_hjs_params * )(p->params);
*)

AA :=  0.757211:
BB := -0.106364:
CC := -0.118649:
DD :=  0.609650:

kF := (rs, z) -> (3.0*Pi^2*(1.0 + z))^(1.0/3.0) * RS_FACTOR/rs:
nu := (rs, z) -> params_a_omega/kF(rs, z):

fH := s -> add(params_a_a[i]*s^(1+i), i=1..6)/(1.0 + add(params_a_b[i]*s^i, i=1..9)):

mzeta   := s -> s^2*fH(s):
meta    := s -> AA + mzeta(s):
mlambda := s -> DD + mzeta(s):
mchi    := (rs, z, s) -> nu(rs, z)/sqrt(mlambda(s) + nu(rs, z)^2):

fF := (rs, z, s) ->
  1.0 - s^2/(27.0*CC*(1.0 + s^2/4.0)) - mzeta(s)/(2.0*CC):

fG := (rs, z, s) ->
  - 2.0/5.0  * CC*fF(rs, z, s)*mlambda(s)
  - 4.0/15.0 * BB*mlambda(s)^2
  - 6.0/5.0  * AA*mlambda(s)^3
  - mlambda(s)^(7.0/2.0)*(4.0/5.0*sqrt(Pi) + 12.0/5.0*(sqrt(mzeta(s)) - sqrt(meta(s)))):

f1 := (rs, z, s) -> 
   + AA 
   - 4.0/9.0 * BB*(1.0 - mchi(rs, z, s))/mlambda(s)
   - 2.0/9.0 * CC*fF(rs, z, s)*(2.0 - 3.0*mchi(rs, z, s) + mchi(rs, z, s)^3)/mlambda(s)^2
   - 1.0/9.0 * fG(rs, z, s)*(8.0 - 15.0*mchi(rs, z, s) + 10.0*mchi(rs, z, s)^3 - 3.0*mchi(rs, z, s)^5)/mlambda(s)^3
   + 2.0*nu(rs, z)*(sqrt(mzeta(s) + nu(rs, z)^2) -  sqrt(meta(s) + nu(rs, z)^2))
   + 2.0*mzeta(s)*log((nu(rs, z) + sqrt(mzeta(s) + nu(rs, z)^2))/(nu(rs, z) + sqrt(mlambda(s) + nu(rs, z)^2)))
   - 2.0*meta(s)*log((nu(rs, z) + sqrt(meta(s) + nu(rs, z)^2))/(nu(rs, z) + sqrt(mlambda(s) + nu(rs, z)^2))):

f_lda := (rs, z) -> -X_FACTOR_C*RS_FACTOR*((1.0 + z)/2.0)^(4.0/3.0)/rs:

f := (rs, z, xt, xs_0_, xs_1_) -> 
  f_lda(rs, z)*f1(rs, z, X2S*xs_0_) + f_lda(rs, -z)*f1(rs, -z, X2S*xs_1_):
