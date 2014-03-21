## Process this file with automake to produce Makefile.in

## Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU Lesser General Public License as published by
## the Free Software Foundation; either version 2, or (at your option)
## any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU Lesser General Public License for more details.
##
## You should have received a copy of the GNU Lesser General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
## 02110-1301, USA.
##
## $Id$

noinst_PROGRAMS = xc-get_data xc-consistency xc-regression xc-error
dist_noinst_SCRIPTS = xc-run_testsuite xc-reference.pl
TESTS = xc-run_testsuite
TESTS_ENVIRONMENT = buildir=$(top_builddir)/testsuite

xc_get_data_SOURCES = xc-get_data.c
xc_get_data_LDADD = -L../src/ -lxc -lm
xc_get_data_CPPFLAGS = -I$(srcdir)/../src/ -I$(top_builddir)/src

xc_consistency_SOURCES = xc-consistency.c
xc_consistency_LDADD = -L../src/ -lxc -lm
xc_consistency_CPPFLAGS = -I$(srcdir)/../src/ -I$(top_builddir)/src

xc_regression_SOURCES = xc-regression.c
xc_regression_LDADD = -L../src/ -lxc -lm
xc_regression_CPPFLAGS = -I$(srcdir)/../src/ -I$(top_builddir)/src

xc_error_SOURCES = xc-error.c
xc_error_LDADD = -lm
xc_error_CPPFLAGS = -I$(srcdir)/../src/ -I$(top_builddir)/src

dist_noinst_DATA = $(srcdir)/df_repo/*.data $(srcdir)/regression/*.pol.bz2 $(srcdir)/regression/*.unpol.bz2 $(srcdir)/input/*

CLEANFILES = *~ *.bak
