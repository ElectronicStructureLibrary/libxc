(* Replace certain mathematica constructions by regular C code *)

cstring = StringReplace[string, 
  {"Power(E," -> "EXP(",
   "Power" -> "POW",
   "Sqrt" -> "SQRT",
   "Pi" -> "M_PI",
   "XFACTORC" -> "X_FACTOR_C",
   "ArcSinh" -> "asinh",
   "Tanh" -> "TANH",
   "PolyLog(2," -> "XC(dilogarithm)("
  }
]
