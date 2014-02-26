/*
 Copyright (C) 2014 Susi Lehtola

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation; either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/


#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include "xc.h"

const char *get_kind(const xc_func_type *func) {
  switch(func->info->kind) {
  case(XC_EXCHANGE):
    return "XC_EXCHANGE";

  case(XC_CORRELATION):
    return "XC_CORRELATION";

  case(XC_EXCHANGE_CORRELATION):
   return "XC_EXCHANGE_CORRELATION";

  case(XC_KINETIC):
   return "XC_KINETIC";

  default:
    printf("Internal error in get_kind.\n");
    return "";
  }
}

const char *get_family(const xc_func_type *func) {
  switch(func->info->family) {
  case(XC_FAMILY_UNKNOWN):
    return "XC_FAMILY_UNKNOWN";

  case(XC_FAMILY_LDA):
    return "XC_FAMILY_LDA";

  case(XC_FAMILY_GGA):
    return "XC_FAMILY_GGA";

  case(XC_FAMILY_MGGA):
    return "XC_FAMILY_MGGA";

  case(XC_FAMILY_LCA):
    return "XC_FAMILY_LCA";

  case(XC_FAMILY_OEP):
    return "XC_FAMILY_OEP";

  case(XC_FAMILY_HYB_GGA):
    return "XC_FAMILY_HYB_GGA";

  case(XC_FAMILY_HYB_MGGA):
    return "XC_FAMILY_HYB_MGGA";

   default:
    printf("Internal error in get_family.\n");
    return "";
  }
}

int main(int argc, char **argv) {
  if(argc!=2) {
    printf("Usage: %s func_id\n",argv[0]);
    return 1;
  }

  int func_id, error;
  xc_func_type func;
  char *fname;

  /* Is functional defined by a string constant? */
  if(isalpha(argv[1][0]))
    func_id=XC(functional_get_number)(argv[1]);
  else
    func_id=atoi(argv[1]);

  /* Initialize functional */
  error=xc_func_init(&func, func_id, XC_UNPOLARIZED);
  if(error) {
    printf("Functional '%s' not found.\n",argv[1]);
    return 1;
  }

  /* Get functional name */
  fname=XC(functional_get_name)(func_id);

  /* Print out info */
  printf("%10s: %-20i\t%10s: %-25s\n","func_id",func_id,"name",fname);
  printf("%10s: %-20s\t%10s: %-25s\n","family",get_family(&func),"kind",get_kind(&func));
  printf("%10s: %s\n","comment",func.info->name);

  printf("\nReference(s):\n");
  printf("%s\n",func.info->refs);

  /* Free memory */
  xc_func_end(&func);
  free(fname);

  return 0;
}
