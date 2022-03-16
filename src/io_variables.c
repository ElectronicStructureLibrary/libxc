/*
 Copyright (C) 2022 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"


/* allocates the output variables on request from the user. The memory
   is initialized to zero. Input parameters are
  np:     Number of grid points
  orders: Array size XC_MAXIMUM_ORDER that determines if orders[i] should be allocated
  family: Family of functional
  flags:  Flags of functional. Necessary to know if LAPL or TAU variables should
          be allocated
  spin:   XC_UNPOLARIZED, XC_POLARIZED
*/

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

xc_output_variables *
xc_output_variables_allocate(double np, const int *orders, int family, int flags, int nspin){
  xc_output_variables *out;
  const xc_dimensions *dim;
  int i;
  
  // initialize the dimension structure
  if(nspin == XC_UNPOLARIZED)
    dim = &dimensions_unpolarized;
  else
    dim = &dimensions_polarized;
  
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
    out->fields[i] = (double *) libxc_malloc(sizeof(double)*np*dim->fields[i+5]);
  }
    
  return out;
}
  
void
xc_output_variables_deallocate(xc_output_variables *out)
{
  int i;
  
  for(i=0; i<XC_TOTAL_NUMBER_OUTPUT_VARIABLES; i++)
    if(out->fields[i] != NULL)
      libxc_free(out->fields[i]);
  libxc_free(out);
}

void
xc_output_variables_initialize(xc_output_variables *out, int np, int nspin)
{
  int i;
  const xc_dimensions *dim;

  // initialize the dimension structure
  if(nspin == XC_UNPOLARIZED)
    dim = &dimensions_unpolarized;
  else
    dim = &dimensions_polarized;
  
  for(i=0; i<XC_TOTAL_NUMBER_OUTPUT_VARIABLES; i++)
    if(out->fields[i] != NULL)
      libxc_memset(out->fields[i], 0, sizeof(double)*np*dim->fields[i+5]);
}
