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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "config.h"

#include "xc.h"

/* xc_config.h needs to be included to use FLOAT and related macros*/
#include "xc_config.h"

int XC(func_info_get_number)(const XC(func_info_type) *info)
{
  return info->number;
}

int XC(func_info_get_kind)(const XC(func_info_type) *info)
{
  return info->kind;
}

char const *XC(func_info_get_name)(const XC(func_info_type) *info)
{
  return info->name;
}

int XC(func_info_get_family)(const XC(func_info_type) *info)
{
  return info->family;
}

int XC(func_info_get_flags)(const XC(func_info_type) *info)
{
  return info->flags;
}

char const *XC(func_info_get_ref)(const XC(func_info_type) *info, int number)
{
  assert(number >=0 && number < 5);

  if (info->refs[number] == NULL) {
    return NULL;
  } else {
    return info->refs[number]->ref;
  }
}
