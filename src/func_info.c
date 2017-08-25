/*
 Copyright (C) 2006-2007 M.A.L. Marques
 Copyright (C) 2016 M. Oliveira

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

#include "util.h"

/* xc_config.h needs to be included to use FLOAT and related macros*/
#include "xc_config.h"

int xc_func_info_get_number(const xc_func_info_type *info)
{
  return info->number;
}

int xc_func_info_get_kind(const xc_func_info_type *info)
{
  return info->kind;
}

char const *xc_func_info_get_name(const xc_func_info_type *info)
{
  return info->name;
}

int xc_func_info_get_family(const xc_func_info_type *info)
{
  return info->family;
}

int xc_func_info_get_flags(const xc_func_info_type *info)
{
  return info->flags;
}

const func_reference_type *xc_func_info_get_references(const xc_func_info_type *info, int number)
{
  assert(number >=0 && number < XC_MAX_REFERENCES);

  if (info->refs[number] == NULL) {
    return NULL;
  } else {
    return info->refs[number];
  }
}

int xc_func_info_get_n_ext_params(xc_func_info_type *info)
{
  assert(info!=NULL);

  return info->n_ext_params;
}

char const *xc_func_info_get_ext_params_description(xc_func_info_type *info, int number)
{
  assert(number >=0 && number < info->n_ext_params);

  return info->ext_params[number].description;
}

double xc_func_info_get_ext_params_default_value(xc_func_info_type *info, int number)
{
  assert(number >=0 && number < info->n_ext_params);

  return info->ext_params[number].value;
}
