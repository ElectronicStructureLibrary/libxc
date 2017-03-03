(* type: work_gga_c *)

$include "gga_c_pw91.mpl"

c1 := 1.1015:
c2 := 0.6625:

f  := (rs, z, xt, xs0, xs1) ->
  + c1*f_pw91(rs, z, xt, xs0, xs1) + (c2 - c1)*(
  +    f_pw91(rs*(2.0*(1.0 + z))^(1.0/3.0),  1.0, xs0, xs0, 0.0)
  +    f_pw91(rs*(2.0*(1.0 - z))^(1.0/3.0), -1.0, xs1, 0.0, xs1)
):