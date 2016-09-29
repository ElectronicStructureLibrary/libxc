(* Replace certain mathematica constructions by regular C code *)

cstring = StringReplace[string, 
  {
    "Pi" -> "M_PI",
    "XFACTORC" -> "X_FACTOR_C",
    "Power(E," -> "EXP(",
    "Power" -> "POW",
    "Sqrt" -> "SQRT",
    "Log"  -> "LOG",
    "ArcSinh" -> "asinh",
    "Tanh" -> "TANH",
    "PolyLog(2," -> "XC(dilogarithm)("
  }
]
