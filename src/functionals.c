#include "xc.h"

extern func_type *lda_known_funct[], *gga_known_funct[];

int family_from_id(int id)
{
  int i;

  /* first let us check if it is an LDA */
  for(i=0; lda_known_funct[i]!=NULL; i++){
    if(lda_known_funct[i]->number == id) return XC_FAMILY_LDA;
  }

  /* or is it a GGA? */
  for(i=0; gga_known_funct[i]!=NULL; i++){
    if(gga_known_funct[i]->number == id) return XC_FAMILY_GGA;
  }

  return XC_FAMILY_UNKNOWN;
}
