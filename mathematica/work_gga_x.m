(* To run the code: math -script "work_gga_x.m" functionals/gga_x_g96.m 3 *)
(* Inputs are: functional; order of the derivative required               *)

(* read functional definition from external file *)
Import[$CommandLine[[4]]]

d[x_] = D[f[x], {x, ToExpression[$CommandLine[[5]]]}]

string = ToString[CForm[Simplify[d[x]]]]

Import["replace.m"]

WriteString["stdout", cstring]
