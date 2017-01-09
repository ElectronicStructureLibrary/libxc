(* type: work_gga_x *)

a5  := 133.983631:
a6  :=   3.217063:
a7  := 136.707378:
a8  :=   3.223476:
a9  :=   2.675484:
a10 :=   3.473804:

f1 := s -> (1.0 - a5*s^a6 + a7*s^a8)/(1 + a9*s^a10):

$include "gga_x_lag.mpl"
unassign('f'):

f := x -> f0(X2S*x) + f1(X2S*x):
