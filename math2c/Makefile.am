## Process this file with automake to produce Makefile.in

## Copyright (C) 2016 M. Marques
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2, or (at your option)
## any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
## 02110-1301, USA.
##
## $Id: Makefile.am 12361 2014-08-06 19:40:53Z dstrubbe $

AM_CPPFLAGS = -I$(top_srcdir) -I$(top_srcdir)/src -I$(top_builddir)/src

noinst_LTLIBRARIES = libmath2c.la

libmath2c_la_SOURCES =  \
        gga_x_g96.c

gga_x_g96.c: xcfuncsh

.PHONY: xcfuncsh
xcfuncsh:
	$(MAKE) -C $(top_builddir)/src xc_funcs.h

