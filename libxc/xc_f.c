#include <stdlib.h>
#include <stdio.h>

#include "xc.h"
#include "config.h"

/* generic interfaces */

void FC_FUNC_(xc_lda_init, XC_LDA_INIT)
		 (void **p, int *functional, int *nspin)
{
	*p = malloc(sizeof(lda_type));
	lda_init((lda_type *)(*p), *functional, *nspin);
}

void FC_FUNC_(xc_lda_end, XC_LDA_END)
		 (void **p)
{
	free(*p);
}

void FC_FUNC_(xc_lda, XC_LDA)
		 (void **p, double *rho, double *e, double *v)
{
	lda((lda_type *)(*p), rho, e, v);
}


/* Now come some special initializations */

/* exchange in the LDA */
void FC_FUNC_(xc_lda_x_init, XC_LDA_X_INIT)
		 (void **p, int *nspin, int *dim, int *rel)
{
	*p = malloc(sizeof(lda_type));
	lda_x_init((lda_type *)(*p), *nspin, *dim, *rel);
}

/* Slater's Xalpha */
void FC_FUNC_(xc_lda_c_xalpha_init, XC_LDA_C_XALPHA_INIT)
		 (void **p, int *nspin, int *dim, int *rel, double *alpha)
{
	*p = malloc(sizeof(lda_type));
	lda_c_xalpha_init((lda_type *)(*p), *nspin, *dim, *rel, *alpha);
}
