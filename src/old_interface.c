/*
 Copyright (C) 2006-2022 M.A.L. Marques
               2021 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

void FUNC( )
  (const xc_func_type *p, size_t np, IN_VARIABLES,
   OUT_VARIABLES_0, OUT_VARIABLES_1, OUT_VARIABLES_2, OUT_VARIABLES_3, OUT_VARIABLES_4)
{
  int order = -1;

  if(zk     != NULL) order = 0;
  if(vrho   != NULL) order = 1;
  if(v2rho2 != NULL) order = 2;
  if(v3rho3 != NULL) order = 3;
  if(v4rho4 != NULL) order = 4;

  if(order < 0) return;

  xc_output_variables out;
  libxc_memset(&out, 0, sizeof(xc_output_variables));

  SET_ORDER_0;
  SET_ORDER_1;
  SET_ORDER_2;
  SET_ORDER_3;
  SET_ORDER_4;

  const xc_input_variables_dimensions *idim = input_variables_dimensions_get(p->nspin);
  const xc_input_variables in = {np, idim, {{INPUT_VARIABLES}}};

  xc_evaluate_func(p, order, &in, &out);
}


/* specializations */
void FUNC(_exc)
  (const xc_func_type *p, size_t np, IN_VARIABLES,
   OUT_VARIABLES_0)
{
  xc_output_variables out;
  libxc_memset(&out, 0, sizeof(xc_output_variables));

  SET_ORDER_0;

  const xc_input_variables_dimensions *idim = input_variables_dimensions_get(p->nspin);
  const xc_input_variables in = {np, idim, {{INPUT_VARIABLES}}};

  xc_evaluate_func(p, 0, &in, &out);
}

void FUNC(_exc_vxc)
  (const xc_func_type *p, size_t np, IN_VARIABLES,
   OUT_VARIABLES_0, OUT_VARIABLES_1)
{
  xc_output_variables out;
  libxc_memset(&out, 0, sizeof(xc_output_variables));

  SET_ORDER_0;
  SET_ORDER_1;

  const xc_input_variables_dimensions *idim = input_variables_dimensions_get(p->nspin);
  const xc_input_variables in = {np, idim, {{INPUT_VARIABLES}}};

  xc_evaluate_func(p, 1, &in, &out);
}

void FUNC(_exc_vxc_fxc)
  (const xc_func_type *p, size_t np, IN_VARIABLES,
    OUT_VARIABLES_0, OUT_VARIABLES_1, OUT_VARIABLES_2)
{
  xc_output_variables out;
  libxc_memset(&out, 0, sizeof(xc_output_variables));

  SET_ORDER_0;
  SET_ORDER_1;
  SET_ORDER_2;

  const xc_input_variables_dimensions *idim = input_variables_dimensions_get(p->nspin);
  const xc_input_variables in = {np, idim, {{INPUT_VARIABLES}}};

  xc_evaluate_func(p, 2, &in, &out);
}

void FUNC(_vxc_fxc)
  (const xc_func_type *p, size_t np, IN_VARIABLES,
   OUT_VARIABLES_1, OUT_VARIABLES_2)
{
  xc_output_variables out;
  libxc_memset(&out, 0, sizeof(xc_output_variables));

  SET_ORDER_1;
  SET_ORDER_2;

  const xc_input_variables_dimensions *idim = input_variables_dimensions_get(p->nspin);
  const xc_input_variables in = {np, idim, {{INPUT_VARIABLES}}};

  xc_evaluate_func(p, 2, &in, &out);
}

void FUNC(_exc_vxc_fxc_kxc)
  (const xc_func_type *p, size_t np, IN_VARIABLES,
   OUT_VARIABLES_0, OUT_VARIABLES_1, OUT_VARIABLES_2, OUT_VARIABLES_3)
{
  xc_output_variables out;
  libxc_memset(&out, 0, sizeof(xc_output_variables));
  
  SET_ORDER_0;
  SET_ORDER_1;
  SET_ORDER_2;
  SET_ORDER_3;

  const xc_input_variables_dimensions *idim = input_variables_dimensions_get(p->nspin);
  const xc_input_variables in = {np, idim, {{INPUT_VARIABLES}}};

  xc_evaluate_func(p, 3, &in, &out);
}

void FUNC(_vxc_fxc_kxc)
  (const xc_func_type *p, size_t np, IN_VARIABLES,
   OUT_VARIABLES_1, OUT_VARIABLES_2, OUT_VARIABLES_3)
{
  xc_output_variables out;
  libxc_memset(&out, 0, sizeof(xc_output_variables));

  SET_ORDER_1;
  SET_ORDER_2;
  SET_ORDER_3;

  const xc_input_variables_dimensions *idim = input_variables_dimensions_get(p->nspin);
  const xc_input_variables in = {np, idim, {{INPUT_VARIABLES}}};

  xc_evaluate_func(p, 3, &in, &out);
}


void FUNC(_vxc)
  (const xc_func_type *p, size_t np, IN_VARIABLES,
   OUT_VARIABLES_1)
{
  xc_output_variables out;
  libxc_memset(&out, 0, sizeof(xc_output_variables));

  SET_ORDER_1;

  const xc_input_variables_dimensions *idim = input_variables_dimensions_get(p->nspin);
  const xc_input_variables in = {np, idim, {{INPUT_VARIABLES}}};

  xc_evaluate_func(p, 1, &in, &out);
}

void FUNC(_fxc)
  (const xc_func_type *p, size_t np, IN_VARIABLES,
   OUT_VARIABLES_2)
{
  xc_output_variables out;
  libxc_memset(&out, 0, sizeof(xc_output_variables));

  SET_ORDER_2;

  const xc_input_variables_dimensions *idim = input_variables_dimensions_get(p->nspin);
  const xc_input_variables in = {np, idim, {{INPUT_VARIABLES}}};

  xc_evaluate_func(p, 2, &in, &out);
}

void FUNC(_mgga_kxc)
  (const xc_func_type *p, size_t np, IN_VARIABLES,
   OUT_VARIABLES_3)
{
  xc_output_variables out;
  libxc_memset(&out, 0, sizeof(xc_output_variables));

  SET_ORDER_3;

  const xc_input_variables_dimensions *idim = input_variables_dimensions_get(p->nspin);
  const xc_input_variables in = {np, idim, {{INPUT_VARIABLES}}};

  xc_evaluate_func(p, 3, &in, &out);
}

void FUNC(_lxc)
  (const xc_func_type *p, size_t np, IN_VARIABLES,
   OUT_VARIABLES_4)
{
  xc_output_variables out;
  libxc_memset(&out, 0, sizeof(xc_output_variables));

  SET_ORDER_4;

  const xc_input_variables_dimensions *idim = input_variables_dimensions_get(p->nspin);
  const xc_input_variables in = {np, idim, {{INPUT_VARIABLES}}};

  xc_evaluate_func(p, 4, &in, &out);
}
