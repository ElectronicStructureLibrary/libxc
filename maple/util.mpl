(* a series of useful definitions *)

RS_FACTOR   := (3.0/(4.0*Pi))^(1.0/3.0):
X2S         := 1.0/(2.0*(6.0*Pi^2)^(1.0/3.0)):
X2S_2D      := 1.0/(2.0*(4.0*Pi)^(1.0/2.0)):

X_FACTOR_C  := 3.0/8.0*(3.0/Pi)^(1.0/3.0)*4.0^(2.0/3.0):

FZETAFACTOR := 2.0^(4.0/3.0) - 2.0:

MU_GE       := 10.0/81.0:
MU_PBE      := 0.06672455060314922*Pi^2/3.0:
KAPPA_PBE   := 0.8040:

(* generic conversion functions *)

n_total    := rs -> (RS_FACTOR/rs)^3:
n_spin     := (rs, z) -> (1.0 + z)*n_total(rs)/2.0:
sigma_spin := (rs, z, xs) -> xs^2*n_spin(rs, z)^(8.0/3.0):
