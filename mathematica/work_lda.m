
(* To run the code: math -script "work_gga_x.m" functionals/gga_x_g96.m 3 *)
(* Inputs are: functional; order of the derivative required               *)

(* read functional definition from external file *)
Import[$CommandLine[[4]]]

If[ToExpression[$CommandLine[[6]]] < 0, 
   {
     flimit[rs_] = Limit[f[rs, zeta], zeta -> 0];
     d[rs_] = D[flimit[rs], {rs, ToExpression[$CommandLine[[5]]]}];
     cstring = ToString[CForm[Simplify[d[varrs]]]];
   },{
    d[rs_, zeta_] = D[D[f[rs, zeta], {rs, ToExpression[$CommandLine[[5]]]}], {zeta, ToExpression[$CommandLine[[6]]]}];
    cstring = ToString[CForm[Simplify[d[varrs, varzeta]]]];
  }
]

WriteString["stdout", cstring]
