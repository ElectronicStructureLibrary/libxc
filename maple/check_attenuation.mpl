(* 2020-12-08 Susi Lehtola

   This script can be used to check the convergence of the series
   expansions of the various attenuation functions.
*)

$include "attenuation.mpl"

(* Ensure we're running in ~ double precision *)
Digits := 16;

check_asymptotics := proc(f, a, expansion_order, padding_order);
  (* Calculate large-a expansion *)
  f_large := a -> eval(throw_out_large_n(convert(series(f(b), b=infinity, expansion_order+padding_order), polynom), expansion_order), b=a):

  (* Store cutoff *)
  readlib(Maple_floats):
  DBL_EPSILON := evalhf(DBL_EPSILON):

  (* The first term is usually inva^2, which gives an expansion of the type
             inva^2(O(1) + O(inva^2) + O(inva^4) + ...)
     meaning that the last term in the bracket is O(inva^{expansion_order-2}).
     This term is numerically zero for a larger than
  *)
  a_cutoff := DBL_EPSILON^(-1/(expansion_order-2)):

  (* Print out the value *)
    printf("Expansion order %3d padding %2d original % e asymptotic % e error % "
           "e\n", expansion_order, padding_order, f(a_cutoff), f_large(a_cutoff), f(a_cutoff)/f_large(a_cutoff)-1);
end proc:

printf("attenuation_erf, pad 10\n");
for eo from 4 by 2 to 100 do
  check_asymptotics(attenuation_erf0, a, eo, 10);
end do;

printf("\nattenuation_erf, pad 20\n");
for eo from 4 by 2 to 100 do
  check_asymptotics(attenuation_erf0, a, eo, 20);
end do;

printf("\nattenuation_yukawa, pad 10\n");
for eo from 4 by 2 to 100 do
  check_asymptotics(attenuation_yukawa0, a, eo, 10);
end do;

printf("\nattenuation_yukawa, pad 20\n");
for eo from 4 by 2 to 100 do
  check_asymptotics(attenuation_yukawa0, a, eo, 20);
end do;

printf("\nattenuation_erf_f2, pad 10\n");
for eo from 4 by 2 to 100 do
  check_asymptotics(attenuation_erf_f20, a, eo, 10);
end do;

printf("\nattenuation_erf_f2, pad 20\n");
for eo from 4 by 2 to 100 do
  check_asymptotics(attenuation_erf_f20, a, eo, 20);
end do;

printf("\nattenuation_erf_f3, pad 10\n");
for eo from 8 by 2 to 100 do
  check_asymptotics(attenuation_erf_f30, a, eo, 10);
end do;

printf("\nattenuation_erf_f3, pad 20\n");
for eo from 8 by 2 to 100 do
  check_asymptotics(attenuation_erf_f30, a, eo, 20);
end do;
