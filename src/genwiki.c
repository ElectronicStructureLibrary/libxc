/*
  Copyright (C) 2015 Susi Lehtola

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
#include <string.h>
#include "xc.h"

#define FUNCMAX 1000

int main(void) {
  int func_id, error, i;
  xc_func_type func;
  char *fname;
  int typeprint;

  /* Families to print out */
  const int families[]={XC_FAMILY_LDA, XC_FAMILY_GGA, XC_FAMILY_HYB_GGA, XC_FAMILY_MGGA, XC_FAMILY_HYB_MGGA};
  const char *famleg[]={"LDA","GGA","hybrid GGA","meta-GGA","hybrid meta-GGA"};
  const int Nfam=sizeof(families)/sizeof(families[0]);
  int ifam;
  /* Types to print out */
  const int types[]={XC_EXCHANGE, XC_CORRELATION, XC_EXCHANGE_CORRELATION, XC_KINETIC};
  const char *typeleg[]={"exchange","correlation","exchange-correlation","kinetic"};
  const int Ntype=sizeof(types)/sizeof(types[0]);
  int itype;

  /* Loop over families, and then types */
  for(ifam=0;ifam<Nfam;ifam++) {
    /* Print wiki section header */
    printf("\n== %s functionals ==\n",famleg[ifam]);
    for(itype=0;itype<Ntype;itype++) {
      typeprint=0;
      
      /* Loop over functional ids */
      for(func_id=0;func_id<FUNCMAX;func_id++) {
	/* Initialize functional */
	error = xc_func_init(&func, func_id, XC_UNPOLARIZED);
	if(error) {
	  /* Functional does not exist */
	  // printf("Functional %i doesn't exist\n",i);
	  continue;
	}
	
	/* Check family and type */
	if(func.info->family != families[ifam] || func.info->kind != types[itype]) {
	  //printf("Functional %i type doesn't match\n",i);
	  XC(func_end)(&func);
	  continue;
	}

	/* Get functional keyword */
	fname = XC(functional_get_name)(func_id);
	/* Convert functional keyword to upper case */
	for(i=0; i<strlen(fname); i++)
	  fname[i]=toupper(fname[i]);
	
	/* Print wiki section header */
	if(!typeprint) {
	  typeprint=1;
	  printf("\n===== %s %s =====\n",famleg[ifam],typeleg[itype]);
	}
	
	/* Print out info */
	printf("* '''%s''': %s\n",fname,func.info->name);
	for(i=0; i<5; i++){
	  if(func.info->refs[i]==NULL) break;
	  if(strlen(func.info->refs[i]->doi) > 0)
	    printf("** [http://dx.doi.org/%s %s] (doi: %s)\n", func.info->refs[i]->doi, func.info->refs[i]->ref, func.info->refs[i]->doi);
	  else
	    printf("** %s\n", func.info->refs[i]->ref);
	}
	
	/* Free memory */
	XC(func_end)(&func);
	free(fname);
      }
    }
  }
  
  return 0;
}
