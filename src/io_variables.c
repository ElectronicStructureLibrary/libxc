/*
 Copyright (C) 2022 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

/* Output variables */

/* mapping input variable -> name */
const char *xc_input_variables_name[XC_TOTAL_NUMBER_INPUT_VARIABLES] =
  {"rho", "sigma", "lapl", "tau", "exx"};

/* mapping input variable -> family */
const int xc_input_variables_family_key[XC_TOTAL_NUMBER_INPUT_VARIABLES] =
  {XC_FAMILY_LDA, XC_FAMILY_GGA, XC_FAMILY_MGGA, XC_FAMILY_MGGA, XC_FAMILY_HGGA};

/* mapping input variable -> flags */
const int xc_input_variables_flags_key[XC_TOTAL_NUMBER_INPUT_VARIABLES] =
  {0, 0, XC_FLAGS_NEEDS_LAPLACIAN, XC_FLAGS_NEEDS_TAU, 0};

const xc_input_variables_dimensions input_variables_dimensions_unpolarized = 
  {1, 1, 1, 1, 1};
const xc_input_variables_dimensions input_variables_dimensions_polarized = 
  {2, 3, 2, 2, 2};

const xc_input_variables_dimensions *input_variables_dimensions_get(int nspin)
{
  if(nspin == XC_UNPOLARIZED)
    return &input_variables_dimensions_unpolarized;
  return &input_variables_dimensions_polarized;
}

/* allocates the input variables on request from the user. The memory
   is not initialized. Input parameters are
  np:     Number of grid points
  family: Family of functional
  flags:  Flags of functional. Necessary to know if LAPL or TAU variables should
          be allocated
  spin:   XC_UNPOLARIZED, XC_POLARIZED
*/
xc_input_variables *
xc_input_variables_allocate(double np, int family, int flags, int nspin){
  xc_input_variables *in;
  int ii;

  /* allocate output structure */
  in = (xc_input_variables *)libxc_malloc(sizeof(xc_input_variables));

  /* initialize the structure to NULLs */
  libxc_memset(in, 0, sizeof(xc_input_variables));

  /* determined spin dimensions */
  in->dim = input_variables_dimensions_get(nspin);

  /* if np == 0 then do not allocate the internal pointers */
  if (np <= 0)
    return in;
  in->np = np;

  for(ii=0; ii<XC_TOTAL_NUMBER_INPUT_VARIABLES; ii++){
    if(family < xc_input_variables_family_key[ii])
      continue;

    if(family >= XC_FAMILY_MGGA){
      if(! (flags & XC_FLAGS_NEEDS_LAPLACIAN) &&
         (xc_input_variables_flags_key[ii] & XC_FLAGS_NEEDS_LAPLACIAN))
        continue;
      if(! (flags & XC_FLAGS_NEEDS_TAU) &&
         (xc_input_variables_flags_key[ii] & XC_FLAGS_NEEDS_TAU))
        continue;
    }
    in->fields[ii] = (double *) libxc_malloc(sizeof(double)*np*in->dim->fields[ii]);
  }
    
  return in;
}

/* 
  This function returns -1 if all relevant variables (as determined by
  orders, family, and flags) are allocated, or the number of the first
  unallocated field otherwise.
*/
int
xc_input_variables_sanity_check(const xc_input_variables *in, int family, int flags)
{
  int ii;

  for(ii=0; ii<XC_TOTAL_NUMBER_INPUT_VARIABLES; ii++){
    if(family < xc_input_variables_family_key[ii])
      continue;

    if(family >= XC_FAMILY_MGGA){
      if(! (flags & XC_FLAGS_NEEDS_LAPLACIAN) &&
         (xc_input_variables_flags_key[ii] & XC_FLAGS_NEEDS_LAPLACIAN))
        continue;
      if(! (flags & XC_FLAGS_NEEDS_TAU) &&
         (xc_input_variables_flags_key[ii] & XC_FLAGS_NEEDS_TAU))
        continue;
    }
    if(in->fields[ii] == NULL)
      return ii; /* the field ii is not allocated */
  }

  /* all relevant fields are allocated */
  return -1;
}


void
xc_input_variables_deallocate(xc_input_variables *in)
{
  int ii;
  
  for(ii=0; ii<XC_TOTAL_NUMBER_INPUT_VARIABLES; ii++)
    if(in->fields[ii] != NULL)
      libxc_free(in->fields[ii]);
  libxc_free(in);
}


void
xc_input_variables_initialize(xc_input_variables *in)
{
  int ii;

  for(ii=0; ii<XC_TOTAL_NUMBER_INPUT_VARIABLES; ii++)
    if(in->fields[ii] != NULL)
      libxc_memset(in->fields[ii], 0, sizeof(double)*in->np*in->dim->fields[ii]);
}


/* Output variables */

/* mapping output variable -> name */
const char *xc_output_variables_name[XC_TOTAL_NUMBER_OUTPUT_VARIABLES] =
  {
   /* order 0 (1 var) */
   "zk",    
   /* order 1 (5 vars) */
   "vrho", "vsigma", "vlapl", "vtau", "vexx",    
   /* order 2 (15 vars) */
   "v2rho2", "v2rhosigma", "v2rholapl", "v2rhotau", "v2rhoexx",
   "v2sigma2", "v2sigmalapl", "v2sigmatau", "v2sigmaexx",
   "v2lapl2", "v2lapltau", "v2laplexx",
   "v2tau2", "v2tauexx",
   "v2exx2",
   /* order 3 (35 vars) */
   "v3rho3", "v3rho2sigma", "v3rho2lapl", "v3rho2tau", "v3rho2exx",
   "v3rhosigma2", "v3rhosigmalapl", "v3rhosigmatau", "v3rhosigmaexx",
   "v3rholapl2", "v3rholapltau", "v3rholaplexx",
   "v3rhotau2", "v3rhotauexx",
   "v3rhoexx2",
   "v3sigma3", "v3sigma2lapl", "v3sigma2tau", "v3sigma2exx",
   "v3sigmalapl2", "v3sigmalapltau", "v3sigmalaplexx",
   "v3sigmatau2", "v3sigmatauexx",
   "v3sigmaexx2",
   "v3lapl3", "v3lapl2tau", "v3lapl2exx",
   "v3lapltau2", "v3lapltauexx",
   "v3laplexx2",
   "v3tau3", "v3tau2exx", "v3tauexx2", "v3exx3", 
   /* order 4 (68 vars) */
   "v4rho4", "v4rho3sigma", "v4rho3lapl", "v4rho3tau", "v4rho3exx",
   "v4rho2sigma2", "v4rho2sigmalapl", "v4rho2sigmatau", "v4rho2sigmaexx",
   "v4rho2lapl2", "v4rho2lapltau", "v4rho2laplexx",
   "v4rho2tau2", "v4rho2tauexx",
   "v4rho2exx2",
   "v4rhosigma3", "v4rhosigma2lapl", "v4rhosigma2tau", "tv4rhosigma2exx",
   "v4rhosigmalapl2", "v4rhosigmalapltau", "v4rhosigmalaplexx",
   "v4rhosigmatau2", "v4rhosigmatauexx",
   "v4rhosigmaexx2",
   "v4rholapl3", "v4rholapl2tau", "v4rholapl2exx",
   "v4rholapltau2", "v4rholapltauexx",
   "v4rholaplexx2",
   "v4rhotau3", "v4rhotau2exx", "v4rhoexx3",
   "v4sigma4", "v4sigma3lapl", "v4sigma3tau", "v4sigma3exx",
   "v4sigma2lapl2", "v4sigma2lapltau", "v4sigma2laplexx",
   "v4sigma2tau2", "v4sigma2tauexx",
   "v4sigma2exx2",
   "v4sigmalapl3", "v4sigmalapl2tau", "v4sigmalapl2exx",
   "v4sigmalapltau2", "v4sigmalapltauexx",
   "v4sigmalaplexx2",
   "v4sigmatau3", "v4sigmatau2exx", "v4sigmatauexx2", "v4sigmaexx3",
   "v4lapl4", "v4lapl3tau", "v4lapl3exx",
   "v4lapl2tau2", "v4lapl2tauexx", "v4lapl2exx2",
   "v4lapltau3", "v4lapltau2exx", "v4lapltauexx2", "v4laplexx3",
   "v4tau4", "v4tau3exx", "v4tauexx3", "v4exx4",
  };

/* mapping output variable -> order of derivative */
const int xc_output_variables_order_key[XC_TOTAL_NUMBER_OUTPUT_VARIABLES] =
  {/* order 0 (1 var) */
   0,
   /* order 1 (5 vars) */
   1, 1, 1, 1, 1,
   /* order 2 (15 vars) */
   2, 2, 2, 2, 2,
   2, 2, 2, 2,
   2, 2, 2,
   2, 2,
   2,
   /* order 3 (35 vars) */
   3, 3, 3, 3, 3,
   3, 3, 3, 3,
   3, 3, 3,
   3, 3,
   3,
   3, 3, 3, 3,
   3, 3, 3,
   3, 3,
   3,
   3, 3, 3,
   3, 3,
   3,
   3, 3, 3, 3,
   /* order 4 (68 vars) */
   4, 4, 4, 4, 4,
   4, 4, 4, 4,
   4, 4, 4,
   4, 4,
   4,
   4, 4, 4, 4,
   4, 4, 4,
   4, 4,
   4,
   4, 4, 4,
   4, 4,
   4,
   4, 4, 4,
   4, 4, 4, 4,
   4, 4, 4,
   4, 4,
   4,
   4, 4, 4,
   4, 4,
   4,
   4, 4, 4, 4,
   4, 4, 4,
   4, 4, 4,
   4, 4, 4, 4,
   4, 4, 4, 4
  };

/* mapping output variable -> family */
const int xc_output_variables_family_key[XC_TOTAL_NUMBER_OUTPUT_VARIABLES] =
  {/* order 0 (1 var) */
   XC_FAMILY_LDA,
   /* order 1 (5 vars) */
   XC_FAMILY_LDA, XC_FAMILY_GGA, XC_FAMILY_MGGA, XC_FAMILY_MGGA, XC_FAMILY_HGGA,
   /* order 2 (15 vars) */
   XC_FAMILY_LDA, XC_FAMILY_GGA, XC_FAMILY_MGGA, XC_FAMILY_MGGA, XC_FAMILY_HGGA,
   XC_FAMILY_GGA, XC_FAMILY_MGGA, XC_FAMILY_MGGA, XC_FAMILY_HGGA,
   XC_FAMILY_MGGA, XC_FAMILY_MGGA, XC_FAMILY_HGGA,
   XC_FAMILY_MGGA, XC_FAMILY_HGGA,
   XC_FAMILY_HGGA,
   /* order 3 (35 vars) */
   XC_FAMILY_LDA, XC_FAMILY_GGA, XC_FAMILY_MGGA, XC_FAMILY_MGGA, XC_FAMILY_HGGA,
   XC_FAMILY_GGA, XC_FAMILY_MGGA, XC_FAMILY_MGGA, XC_FAMILY_HGGA,
   XC_FAMILY_MGGA, XC_FAMILY_MGGA, XC_FAMILY_HGGA,
   XC_FAMILY_MGGA, XC_FAMILY_HGGA,
   XC_FAMILY_HGGA,
   XC_FAMILY_GGA, XC_FAMILY_MGGA, XC_FAMILY_MGGA, XC_FAMILY_HGGA,
   XC_FAMILY_MGGA, XC_FAMILY_MGGA, XC_FAMILY_HGGA,
   XC_FAMILY_MGGA, XC_FAMILY_HGGA,
   XC_FAMILY_HGGA,
   XC_FAMILY_MGGA, XC_FAMILY_MGGA, XC_FAMILY_HGGA,
   XC_FAMILY_MGGA, XC_FAMILY_HGGA,
   XC_FAMILY_HGGA,
   XC_FAMILY_MGGA, XC_FAMILY_HGGA, XC_FAMILY_HGGA, XC_FAMILY_HGGA,
   /* order 4 (68 vars) */
   XC_FAMILY_LDA, XC_FAMILY_GGA, XC_FAMILY_MGGA, XC_FAMILY_MGGA, XC_FAMILY_HGGA,
   XC_FAMILY_GGA, XC_FAMILY_MGGA, XC_FAMILY_MGGA, XC_FAMILY_HGGA,
   XC_FAMILY_MGGA, XC_FAMILY_MGGA, XC_FAMILY_HGGA,
   XC_FAMILY_MGGA, XC_FAMILY_HGGA,
   XC_FAMILY_HGGA,
   XC_FAMILY_GGA, XC_FAMILY_MGGA, XC_FAMILY_MGGA, XC_FAMILY_HGGA,
   XC_FAMILY_MGGA, XC_FAMILY_MGGA, XC_FAMILY_HGGA,
   XC_FAMILY_MGGA, XC_FAMILY_HGGA,
   XC_FAMILY_HGGA,
   XC_FAMILY_MGGA, XC_FAMILY_MGGA, XC_FAMILY_HGGA,
   XC_FAMILY_MGGA, XC_FAMILY_HGGA,
   XC_FAMILY_HGGA,
   XC_FAMILY_MGGA, XC_FAMILY_HGGA, XC_FAMILY_HGGA,
   XC_FAMILY_GGA, XC_FAMILY_MGGA, XC_FAMILY_MGGA, XC_FAMILY_HGGA,
   XC_FAMILY_MGGA, XC_FAMILY_MGGA, XC_FAMILY_HGGA,
   XC_FAMILY_MGGA, XC_FAMILY_HGGA,
   XC_FAMILY_HGGA,
   XC_FAMILY_MGGA, XC_FAMILY_MGGA, XC_FAMILY_HGGA,
   XC_FAMILY_MGGA, XC_FAMILY_HGGA,
   XC_FAMILY_HGGA,
   XC_FAMILY_MGGA, XC_FAMILY_HGGA, XC_FAMILY_HGGA, XC_FAMILY_HGGA,
   XC_FAMILY_MGGA, XC_FAMILY_MGGA, XC_FAMILY_HGGA,
   XC_FAMILY_MGGA, XC_FAMILY_HGGA, XC_FAMILY_HGGA,
   XC_FAMILY_MGGA, XC_FAMILY_HGGA, XC_FAMILY_HGGA, XC_FAMILY_HGGA,
   XC_FAMILY_MGGA, XC_FAMILY_HGGA, XC_FAMILY_HGGA, XC_FAMILY_HGGA
  };

/* mapping output variable -> flags */
const int xc_output_variables_flags_key[XC_TOTAL_NUMBER_OUTPUT_VARIABLES] =
  {/* order 0 (1 var) */
   0,
   /* order 1 (5 vars) */
   0, 0, XC_FLAGS_NEEDS_LAPLACIAN, XC_FLAGS_NEEDS_TAU, 0,
   /* order 2 (15 vars) */
   0, 0, XC_FLAGS_NEEDS_LAPLACIAN, XC_FLAGS_NEEDS_TAU, 0,
   0, XC_FLAGS_NEEDS_LAPLACIAN, XC_FLAGS_NEEDS_TAU, 0,
   XC_FLAGS_NEEDS_LAPLACIAN, XC_FLAGS_NEEDS_LAPLACIAN | XC_FLAGS_NEEDS_TAU, XC_FLAGS_NEEDS_LAPLACIAN,
   XC_FLAGS_NEEDS_TAU, XC_FLAGS_NEEDS_TAU,
   0,
   /* order 3 (35 vars) */
   0, 0, XC_FLAGS_NEEDS_LAPLACIAN, XC_FLAGS_NEEDS_TAU, 0,
   0, XC_FLAGS_NEEDS_LAPLACIAN, XC_FLAGS_NEEDS_TAU, 0,
   XC_FLAGS_NEEDS_LAPLACIAN, XC_FLAGS_NEEDS_LAPLACIAN | XC_FLAGS_NEEDS_TAU, XC_FLAGS_NEEDS_LAPLACIAN,
   XC_FLAGS_NEEDS_TAU, XC_FLAGS_NEEDS_TAU,
   0,
   0, XC_FLAGS_NEEDS_LAPLACIAN, XC_FLAGS_NEEDS_TAU, 0,
   XC_FLAGS_NEEDS_LAPLACIAN, XC_FLAGS_NEEDS_LAPLACIAN | XC_FLAGS_NEEDS_TAU, XC_FLAGS_NEEDS_LAPLACIAN,
   XC_FLAGS_NEEDS_TAU, XC_FLAGS_NEEDS_TAU,
   0,
   XC_FLAGS_NEEDS_LAPLACIAN, XC_FLAGS_NEEDS_LAPLACIAN | XC_FLAGS_NEEDS_TAU, XC_FLAGS_NEEDS_LAPLACIAN,
   XC_FLAGS_NEEDS_LAPLACIAN | XC_FLAGS_NEEDS_TAU, XC_FLAGS_NEEDS_LAPLACIAN | XC_FLAGS_NEEDS_TAU,
   XC_FLAGS_NEEDS_LAPLACIAN,
   XC_FLAGS_NEEDS_TAU, XC_FLAGS_NEEDS_TAU, XC_FLAGS_NEEDS_TAU, 0,
   /* order 4 (68 vars) */
   0, 0, XC_FLAGS_NEEDS_LAPLACIAN, XC_FLAGS_NEEDS_TAU, 0,
   0, XC_FLAGS_NEEDS_LAPLACIAN, XC_FLAGS_NEEDS_TAU, 0,
   XC_FLAGS_NEEDS_LAPLACIAN, XC_FLAGS_NEEDS_LAPLACIAN | XC_FLAGS_NEEDS_TAU, XC_FLAGS_NEEDS_LAPLACIAN,
   XC_FLAGS_NEEDS_TAU, XC_FLAGS_NEEDS_TAU,
   0,
   0, XC_FLAGS_NEEDS_LAPLACIAN, XC_FLAGS_NEEDS_TAU, 0,
   XC_FLAGS_NEEDS_LAPLACIAN, XC_FLAGS_NEEDS_LAPLACIAN | XC_FLAGS_NEEDS_TAU, XC_FLAGS_NEEDS_LAPLACIAN,
   XC_FLAGS_NEEDS_TAU, XC_FLAGS_NEEDS_TAU,
   0,
   XC_FLAGS_NEEDS_LAPLACIAN, XC_FLAGS_NEEDS_LAPLACIAN | XC_FLAGS_NEEDS_TAU, XC_FLAGS_NEEDS_LAPLACIAN,
   XC_FLAGS_NEEDS_LAPLACIAN | XC_FLAGS_NEEDS_TAU, XC_FLAGS_NEEDS_LAPLACIAN | XC_FLAGS_NEEDS_TAU,
   XC_FLAGS_NEEDS_LAPLACIAN,
   XC_FLAGS_NEEDS_TAU, XC_FLAGS_NEEDS_TAU, 0,
   0, XC_FLAGS_NEEDS_LAPLACIAN, XC_FLAGS_NEEDS_TAU, 0,
   XC_FLAGS_NEEDS_LAPLACIAN, XC_FLAGS_NEEDS_LAPLACIAN | XC_FLAGS_NEEDS_TAU, XC_FLAGS_NEEDS_LAPLACIAN,
   XC_FLAGS_NEEDS_TAU, XC_FLAGS_NEEDS_TAU,
   0,
   XC_FLAGS_NEEDS_LAPLACIAN, XC_FLAGS_NEEDS_LAPLACIAN | XC_FLAGS_NEEDS_TAU, XC_FLAGS_NEEDS_LAPLACIAN,
   XC_FLAGS_NEEDS_LAPLACIAN | XC_FLAGS_NEEDS_TAU, XC_FLAGS_NEEDS_LAPLACIAN | XC_FLAGS_NEEDS_TAU,
   XC_FLAGS_NEEDS_LAPLACIAN,
   XC_FLAGS_NEEDS_TAU, XC_FLAGS_NEEDS_TAU, XC_FLAGS_NEEDS_TAU, 0,
   XC_FLAGS_NEEDS_LAPLACIAN, XC_FLAGS_NEEDS_LAPLACIAN | XC_FLAGS_NEEDS_TAU, XC_FLAGS_NEEDS_LAPLACIAN,
   XC_FLAGS_NEEDS_LAPLACIAN | XC_FLAGS_NEEDS_TAU, XC_FLAGS_NEEDS_LAPLACIAN | XC_FLAGS_NEEDS_TAU, XC_FLAGS_NEEDS_LAPLACIAN,
   XC_FLAGS_NEEDS_LAPLACIAN | XC_FLAGS_NEEDS_TAU, XC_FLAGS_NEEDS_LAPLACIAN | XC_FLAGS_NEEDS_TAU, XC_FLAGS_NEEDS_LAPLACIAN | XC_FLAGS_NEEDS_TAU, XC_FLAGS_NEEDS_LAPLACIAN,
   XC_FLAGS_NEEDS_TAU, XC_FLAGS_NEEDS_TAU, XC_FLAGS_NEEDS_TAU, 0,
  };


/* these are the dimension of the input and output
   arrays for spin unpolarized and polarized */
const xc_output_variables_dimensions output_variables_dimensions_unpolarized =
{
  /* order 0 */
  1,                /* zk */
  /* order 1 */
  1, 1, 1, 1, 1,    /* vrho, vsigma, vlapl, vtau, vexx */
  /* order 2 */
  1, 1, 1, 1, 1,    /* v2rho2, v2rhosigma, v2rholapl, v2rhotau, v2rhoexx */
  1, 1, 1, 1,       /* v2sigma2, v2sigmalapl, v2sigmatau, v2sigmaexx */
  1, 1, 1,          /* v2lapl2, v2lapltau, v2laplexx */
  1, 1,             /* v2tau2, v2tauexx */
  1,                /* v2exx2 */
  /* order 3 */
  1, 1, 1, 1, 1,    /* v3rho3, v3rho2sigma, v3rho2lapl, v3rho2tau, v3rho2exx */
  1, 1, 1, 1,       /* v3rhosigma2, v3rhosigmalapl, v3rhosigmatau, v3rhosigmaexx */
  1, 1, 1,          /* v3rholapl2, v3rholapltau, v3rholaplexx */
  1, 1,             /* v3rhotau2, v3rhotauexx */
  1,                /* v3rhoexx2 */
  1, 1, 1, 1,       /* v3sigma3, v3sigma2lapl, v3sigma2tau, v3sigma2exx */
  1, 1, 1,          /* v3sigmalapl2, v3sigmalapltau, v3sigmalaplexx */
  1, 1,             /* v3sigmatau2, v3sigmatauexx */
  1,                /* v3sigmaexx2 */
  1, 1, 1,          /* v3lapl3, v3lapl2tau, v3lapl2exx */
  1, 1,             /* v3lapltau2, v3lapltauexx */
  1,                /* v3laplexx2 */
  1, 1, 1, 1,       /* v3tau3, v3tau2exx, v3tauexx2, v3exx3 */
  /* order 4 */
  1, 1, 1, 1, 1,    /* v4rho4, v4rho3sigma, v4rho3lapl, v4rho3tau, v4rho3exx */
  1, 1, 1, 1,       /* v4rho2sigma2, v4rho2sigmalapl, v4rho2sigmatau, v4rho2sigmaexx */
  1, 1, 1,          /* v4rho2lapl2, v4rho2lapltau, v4rho2laplexx */
  1, 1,             /* v4rho2tau2, v4rho2tauexx */
  1,                /* v4rho2exx2 */
  1, 1, 1, 1,       /* v4rhosigma3, v4rhosigma2lapl, v4rhosigma2tau, v4rhosigma2exx */
  1, 1, 1,          /* v4rhosigmalapl2, v4rhosigmalapltau, v4rhosigmalaplexx */
  1, 1,             /* v4rhosigmatau2, v4rhosigmatauexx */
  1,                /* v4rhosigmaexx2 */
  1, 1, 1,          /* v4rhola1pl3, v4rholapl2tau, v4rholapl2exx */
  1, 1,             /* v4rholapltau2, v4rholapltauexx */
  1,                /* v4rholaplexx2 */
  1, 1, 1,          /* v4rhotau3, v4rhotau2exx, v4rhoexx3 */
  1, 1, 1, 1,       /* v4sigma4, v4sigma3lapl, v4sigma3tau, v4sigma3exx */
  1, 1, 1,          /* v4sigma2lapl2, v4sigma2lapltau, v4sigma2laplexx */
  1, 1,             /* v4sigma2tau2, v4sigma2tauexx */
  1,                /* v4sigma2exx2 */
  1, 1, 1,          /* v4sigmalapl3, v4sigmalapl2tau, v4sigmalapl2exx */
  1, 1,             /* v4sigmalapltau2, v4sigmalapltauexx */
  1,                /* v4sigmalaplexx2 */
  1, 1, 1, 1,       /* v4sigmatau3, v4sigmatau2exx, v4sigmatauexx2, v4sigmaexx3 */
  1, 1, 1,          /* v4lapl4, v4lapl3tau, v4lapl3exx */
  1, 1, 1,          /* v4lapl2tau2, v4lapl2tauexx, v4lapl2exx2 */
  1, 1, 1, 1,       /* v4lapltau3, v4lapltau2exx, v4lapltauexx2, v4laplexx3 */
  1, 1, 1, 1        /* v4tau4, v4tau3exx, v4tauexx3, v4exx4 */
};

const xc_output_variables_dimensions output_variables_dimensions_polarized =
{
  /* order 0 */
  1,                /* zk */
  /* order 1 */
  2, 3, 2, 2, 2,    /* vrho, vsigma, vlapl, vtau, vexx */
  /* order 2 */
  3, 6, 4, 4, 4,    /* v2rho2, v2rhosigma, v2rholapl, v2rhotau, v2rhoexx */
  6, 6, 6, 6,       /* v2sigma2, v2sigmalapl, v2sigmatau, v2sigmaexx */
  3, 4, 4,          /* v2lapl2, v2lapltau, v2laplexx */
  3, 4,             /* v2tau2, v2tauexx */
  3,                /* v2exx2 */
  /* order 3 */
  4, 9, 6, 6, 6,    /* v3rho3, v3rho2sigma, v3rho2lapl, v3rho2tau, v3rho2exx */
  12, 12, 12, 12,   /* v3rhosigma2, v3rhosigmalapl, v3rhosigmatau, v3rhosigmaexx */
  6, 8, 8,          /* v3rholapl2, v3rholapltau, v3rholaplexx */
  6, 8,             /* v3rhotau2, v3rhotauexx */
  6,                /* v3rhoexx2 */
  10, 12, 12, 12,   /* v3sigma3, v3sigma2lapl, v3sigma2tau, v3sigma2exx */
  9, 12, 12,        /* v3sigmalapl2, v3sigmalapltau, v3sigmalaplexx */
  9, 12,            /* v3sigmatau2, v3sigmatauexx */
  9,                /* v3sigmaexx2 */
  4, 6, 6,          /* v3lapl3, v3lapl2tau, v3lapl2exx */
  6, 8,             /* v3lapltau2, v3lapltauexx */
  6,                /* v3laplexx2 */
  4, 6, 6, 4,       /* v3tau3, v3tau2exx, v3tauexx2, v3exx3 */
  /* order 4 */
  5, 12, 8, 8, 8,   /* v4rho4, v4rho3sigma, v4rho3lapl, v4rho3tau, v4rho3exx */
  18, 18, 18, 18,   /* v4rho2sigma2, v4rho2sigmalapl, v4rho2sigmatau, v4rho2sigmaexx */
  9, 12, 12,        /* v4rho2lapl2, v4rho2lapltau, v4rho2laplexx */
  9, 12,            /* v4rho2tau2, v4rho2tauexx */
  9,                /* v4rho2exx2 */
  20, 36, 36, 36,   /* v4rhosigma3, v4rhosigma2lapl, v4rhosigma2tau, v4rhosigma2exx */
  18, 24, 24,       /* v4rhosigmalapl2, v4rhosigmalapltau, v4rhosigmalaplexx */
  18, 24,           /* v4rhosigmatau2, v4rhosigmatauexx */
  18,               /* v4rhosigmaexx2 */
  8, 12, 12,        /* v4rholapl3, v4rholapl2tau, v4rholapl2exx */
  12, 16,           /* v4rholapltau2, v4rholapltauexx */
  12,               /* v4rholaplexx2 */
  8, 12, 8,         /* v4rhotau3, v4rhotau2exx, v4rhoexx3 */
  15, 20, 20, 20,   /* v4sigma4, v4sigma3lapl, v4sigma3tau, v4sigma3exx */
  18, 24, 24,       /* v4sigma2lapl2, v4sigma2lapltau, v4sigma2laplexx */
  18, 24,           /* v4sigma2tau2, v4sigma2tauexx */
  18,               /* v4sigma2exx2 */
  12, 18, 18,       /* v4sigmalapl3, v4sigmalapl2tau, v4sigmalapl2exx */
  18, 24,           /* v4sigmalapltau2, v4sigmalapltauexx */
  18,               /* v4sigmalaplexx2 */
  12, 18, 18, 12,   /* v4sigmatau3, v4sigmatau2exx, v4sigmatauexx2, v4sigmaexx3 */
  5, 8, 8,          /* v4lapl4, v4lapl3tau, v4lapl3exx */
  9, 12, 9,         /* v4lapl2tau2, v4lapl2tauexx, v4lapl2exx2 */
  8, 12, 12, 8,     /* v4lapltau3, v4lapltau2exx, v4lapltauexx2, v4laplexx3 */
  5, 8, 8, 5        /* v4tau4, v4tau3exx, v4tauexx3, v4exx4 */
};

const xc_output_variables_dimensions *output_variables_dimensions_get(int nspin)
{
  if(nspin == XC_UNPOLARIZED)
    return &output_variables_dimensions_unpolarized;
  return &output_variables_dimensions_polarized;
}

/* allocates the output variables on request from the user. The memory
   is not initialized. Input parameters are
  np:     Number of grid points
  orders: Array size XC_MAXIMUM_ORDER that determines if orders[i] should be allocated
  family: Family of functional
  flags:  Flags of functional. Necessary to know if LAPL or TAU variables should
          be allocated
  spin:   XC_UNPOLARIZED, XC_POLARIZED
*/
xc_output_variables *
xc_output_variables_allocate(double np, const int *orders, int family, int flags, int nspin){
  xc_output_variables *out;
  const xc_output_variables_dimensions *dim = output_variables_dimensions_get(nspin);
  int i;
  
  /* allocate output structure */
  out = (xc_output_variables *)libxc_malloc(sizeof(xc_output_variables));

  /* initialize the structure to NULLs */
  libxc_memset(out, 0, sizeof(xc_output_variables));

  /* if np == 0 then do not allocate the internal pointers */
  if (np <= 0)
    return out;

  for(i=0; i<XC_TOTAL_NUMBER_OUTPUT_VARIABLES; i++){
    if(! orders[xc_output_variables_order_key[i]])
      continue;

    if(family < xc_output_variables_family_key[i])
      continue;

    if(family >= XC_FAMILY_MGGA){
      if(! (flags & XC_FLAGS_NEEDS_LAPLACIAN) &&
         (xc_output_variables_flags_key[i] & XC_FLAGS_NEEDS_LAPLACIAN))
        continue;
      if(! (flags & XC_FLAGS_NEEDS_TAU) &&
         (xc_output_variables_flags_key[i] & XC_FLAGS_NEEDS_TAU))
        continue;
    }
    out->fields[i] = (double *) libxc_malloc(sizeof(double)*np*dim->fields[i]);
  }
    
  return out;
}
  
/* 
  This function returns -1 if all relevant variables (as determined by
  orders, family, and flags) are allocated, or the number of the first
  unallocated field otherwise. It also returns 10000+order if a
  variable was requested that is not available.
*/
int
xc_output_variables_sanity_check(const xc_output_variables *out, const int *orders, int family, int flags)
{
  int ii;

  if(orders[0] && !(flags & XC_FLAGS_HAVE_EXC))
    return 10000;

  if(orders[1] && !(flags & XC_FLAGS_HAVE_VXC))
    return 10001;

  if(orders[2] && !(flags & XC_FLAGS_HAVE_FXC))
    return 10002;

  if(orders[3] && !(flags & XC_FLAGS_HAVE_KXC))
    return 10003;

  if(orders[4] && !(flags & XC_FLAGS_HAVE_LXC))
    return 10004;
  
  for(ii=0; ii<XC_TOTAL_NUMBER_OUTPUT_VARIABLES; ii++){
    if(! orders[xc_output_variables_order_key[ii]])
      continue;

    if(family < xc_output_variables_family_key[ii])
      continue;

    if(family >= XC_FAMILY_MGGA){
      if(! (flags & XC_FLAGS_NEEDS_LAPLACIAN) &&
         (xc_output_variables_flags_key[ii] & XC_FLAGS_NEEDS_LAPLACIAN))
        continue;
      if(! (flags & XC_FLAGS_NEEDS_TAU) &&
         (xc_output_variables_flags_key[ii] & XC_FLAGS_NEEDS_TAU))
        continue;
    }
    if(out->fields[ii] == NULL)
      return ii; /* the field ii is not allocated */
  }

  /* all relevant fields are allocated */
  return -1;
}


void
xc_output_variables_deallocate(xc_output_variables *out)
{
  int ii;
  
  for(ii=0; ii<XC_TOTAL_NUMBER_OUTPUT_VARIABLES; ii++)
    if(out->fields[ii] != NULL)
      libxc_free(out->fields[ii]);
  libxc_free(out);
}

void
xc_output_variables_initialize(xc_output_variables *out, int np, int nspin)
{
  int ii;
  const xc_output_variables_dimensions *dim = output_variables_dimensions_get(nspin);
  
  for(ii=0; ii<XC_TOTAL_NUMBER_OUTPUT_VARIABLES; ii++)
    if(out->fields[ii] != NULL)
      libxc_memset(out->fields[ii], 0, sizeof(double)*np*dim->fields[ii]);
}
