# Minimum and maximum functions that I can differentiate
m_min := (x1, x2) -> x1 + (x2 - x1)*Heaviside(x1 - x2):
m_max := (x1, x2) -> x1 + (x2 - x1)*Heaviside(x2 - x1):
m_abs := (x) -> m_max(x, -x):

# a series of useful definitions

M_C         := 137.0359996287515: (* speed of light *)

RS_FACTOR   := (3.0/(4.0*Pi))^(1.0/3.0):
X2S         := 1.0/(2.0*(6.0*Pi^2)^(1.0/3.0)):
X2S_2D      := 1.0/(2.0*(4.0*Pi)^(1.0/2.0)):

X_FACTOR_C    := 3.0/8.0*(3.0/Pi)^(1.0/3.0)*4.0^(2.0/3.0):
X_FACTOR_2D_C := 8.0/(3.0*sqrt(Pi)):
K_FACTOR_C    := 3.0/10.0*(6.0*Pi^2)^(2.0/3.0):

MU_GE       := 10.0/81.0:
MU_PBE      := 0.06672455060314922*Pi^2/3.0:
KAPPA_PBE   := 0.8040:

# generic conversion functions
n_total    := rs -> (RS_FACTOR/rs)^3:
n_spin     := (rs, z) -> (1.0 + z)*n_total(rs)/2.0:
sigma_spin := (rs, z, xs) -> xs^2*n_spin(rs, z)^(8.0/3.0):
t_total    := (z, ts0, ts1) ->
  (ts0*((1.0 + z)/2.0)^(5.0/3.0) + ts1*((1.0 - z)/2.0)^(5.0/3.0)):
u_total    := (z, us0, us1) -> t_total(z, us0, us1):

# useful formulas that enter several functionals follow

# von Weizsäcker term
t_vw := (z, xt, us0, us1) -> (xt^2 - u_total(z, us0, us1))/8.0:

# See Eq. (9) of Perdew1992_13244
f_zeta    := z -> ((1.0 + z)^(4.0/3.0) + (1.0 - z)^(4.0/3.0) - 2.0)/(2.0^(4.0/3.0) - 2.0):
f_zeta_2d := z -> 0.5*((1.0 + z)^1.5 + (1.0 - z)^1.5):

# used in several correlation functionals
mphi := z -> ((1.0 + z)^(2.0/3.0) + (1.0 - z)^(2.0/3.0))/2.0:
tt   := (rs, z, xt) -> xt/(4.0*2.0^(1.0/3.0)*mphi(z)*sqrt(rs)):

# in the paper it is beta_a = 0.066725
beta_Hu_Langreth := rs -> 0.066724550603149220*(1.0 + 0.1*rs)/(1.0 + 0.1778*rs):

# This is the Stoll decomposition in our language
lda_stoll_par  := (lda_func, rs, z, spin) ->
  lda_func(rs*(2.0/(1.0 + z))^(1.0/3.0), spin)*(1.0 + z)/2.0:

lda_stoll_perp := (lda_func, rs, z) ->
  + lda_func(rs, z)
  - lda_stoll_par(lda_func, rs,  z,  1.0)
  - lda_stoll_par(lda_func, rs, -z, -1.0):

gga_stoll_par  := (gga_func, rs, z, xs, spin) ->
 + gga_func(rs*(2.0/(1.0 + z))^(1.0/3.0), spin, xs, xs*(1.0 + spin)/2.0, xs*(1.0 - spin)/2.0)
 * (1.0 + z)/2.0:

# Becke function used in several correlation functionals
b88_R_F := (f_x, rs, z, xs) ->
  1.0/(2.0*X_FACTOR_C*n_spin(rs, z)^(1/3)*f_x(xs)):

b88_zss := (css, f_x, rs, z, xs) ->
  2.0*css*b88_R_F(f_x, rs, z, xs):
b88_zab := (cab, f_x, rs, z, xs0, xs1) ->
  cab*(b88_R_F(f_x, rs, z, xs0) + b88_R_F(f_x, rs, -z, xs1)):

# Power series often used in mggas
mgga_w := t -> 1.0*(K_FACTOR_C - t)/(K_FACTOR_C + t):
mgga_series_w := (a, n, t) -> add(a[i]*mgga_w(t)^(i-1), i=1..n):

# Used in screened functionals
kF := (rs, z) -> (3.0*Pi^2*(1.0 + z))^(1.0/3.0) * RS_FACTOR/rs:
nu := (rs, z) -> params_a_omega/kF(rs, z):
