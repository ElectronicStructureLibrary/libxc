/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <math.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>

#include <xc.h>
#include <util.h>

int main(int argc, char *argv[])
{
  xc_func_type func;
  xc_input_variables *in;
  xc_output_variables *out;
  const xc_output_variables_dimensions *out_dim;
  
  int functional, nspin, ii, is, ninput;
  int orders[XC_MAXIMUM_ORDER+1] = {1, 1, 1, 0, 0};

  if(argc < 4){
    printf("Usage:\n%s funct pol <input>\n", argv[0]);
    return 1;
  }

  /* Is functional defined by a string constant? */
  if(isalpha(argv[1][0]))
    functional = xc_functional_get_number(argv[1]);
  else
    functional = atoi(argv[1]);

  nspin = atoi(argv[2]);
  if(nspin == XC_UNPOLARIZED)
    printf("Unpolarized calculation\n");
  else if(nspin == XC_POLARIZED)
    printf("Polarized calculation\n");
  else {
    printf("Invalid value for pol input.\n");
    exit(1);
  }

  /* Initialize functional */
  if(xc_func_init(&func, functional, nspin) != 0){
    fprintf(stderr, "Functional '%d' not found\n", functional);
    exit(1);
  }

  /* define and read input variables */
  in = xc_input_variables_allocate(1, func.info->family, func.info->flags, nspin);
  
  /* let us check how many input parameters our functional has */
  for(ninput=0, ii=0; ii<XC_TOTAL_NUMBER_INPUT_VARIABLES; ii++)
    if(in->fields[ii] != NULL) ninput += in->dim->fields[ii];

  if(argc != 3 + ninput){
    fprintf(stderr, "Wrong number of parameters. Functional expects:\n  ");
    for(ii=0; ii<XC_TOTAL_NUMBER_INPUT_VARIABLES; ii++){
      if(in->fields[ii] == NULL)
        continue;
      for(is=0; is<in->dim->fields[ii]; is++)
        fprintf(stderr, " %s%d ", xc_input_variables_name[ii], is);
    }
    exit(1);
  }
  
  /* Read data */
  for(ninput=0, ii=0; ii<XC_TOTAL_NUMBER_INPUT_VARIABLES; ii++){
    if(in->fields[ii] == NULL)
      continue;
    for(is=0; is<in->dim->fields[ii]; is++)
      in->fields[ii][is] = atof(argv[3 + ninput++]);
  }
  
  /* allocate buffers */
  out = xc_output_variables_allocate(1, orders,
                                     func.info->family,
                                     func.info->flags,
                                     nspin);
  /* xc_output_variables_initialize(out); */
  
  xc_evaluate_func(&func, 2, in, out);

  /* transform to energy per volume */
  if(nspin == XC_UNPOLARIZED)
    out->zk[0] *= in->rho[0];
  else
    out->zk[0] *= in->rho[0] + in->rho[1];
  
  /* let us print the input */
  for(ii=0; ii<XC_TOTAL_NUMBER_INPUT_VARIABLES; ii++){
    if(in->fields[ii] == NULL)
      continue;
    for(is=0; is<in->dim->fields[ii]; is++)
      printf(" %s%d = %#0.2E", xc_input_variables_name[ii], is, in->fields[ii][is]);
  }
  printf("\n");
  
  /* and now we print the output */
  out_dim = output_variables_dimensions_get(nspin);
  for(ninput=0, ii=0; ii<XC_TOTAL_NUMBER_OUTPUT_VARIABLES; ii++){
    if(out->fields[ii] == NULL)
      continue;
    for(is=0; is<out_dim->fields[ii]; is++)
      printf("%3d: %20s[%2d] = %#19.12E\n", ninput++,
             xc_output_variables_name[ii], is, out->fields[ii][is]);
  }

  /* clean up */
  xc_input_variables_deallocate(in);
  xc_output_variables_deallocate(out);
  xc_func_end(&func);

  return 0;
}
