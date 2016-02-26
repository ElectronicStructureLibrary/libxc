/*
 Copyright (C) 2012 M.A.L. Marques, M. Oliveira

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

#include <stdio.h>
#include <string.h>

#include "xc.h"
#include "config.h"

void XC(version)(int *major, int *minor, int *micro) {
  const char *version_string = PACKAGE_VERSION;

  *major = -1;
  *minor = -1;
  *micro = -1;
  sscanf(version_string,"%d.%d.%d",major,minor,micro);

}

char *XC(version_string)() {
  return strdup(PACKAGE_VERSION);
}
