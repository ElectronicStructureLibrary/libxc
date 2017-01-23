(* type: work_gga_x *)

$define gga_x_rpw86_params
$include "gga_x_pw86.mpl"

malpha := 0.02178:
mbeta  := 1.15:
muLV   := 0.8491/9.0:

f0 := s -> 
   + (1.0 + muLV*s^2)/(1.0 + malpha*s^6) 
   + malpha*s^6*f0_pw86(s)/(mbeta + malpha*s^6):

f  := x -> f0(X2S*x):

