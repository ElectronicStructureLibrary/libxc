#include <stdlib.h>

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
		 (void **p, double *rho, double *ex, double *vx)
{
	lda((lda_type *)(*p), rho, ex, vx);
}


/* Now come some special initializations */

/* exchange in the LDA */
void FC_FUNC_(xc_lda_x_init, XC_LDA_X_INIT)
		 (void **p, int *nspin, int *dim, int *rel)
{
	*p = malloc(sizeof(lda_type));
	lda_x_init((lda_type *)(*p), *nspin, *dim, *rel);
}

