## Process this file with automake to produce Makefile.in

## Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
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
## Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
## 02111-1307, USA.
##
## $Id:  $

check_PROGRAMS = xc-get_data xc-consistency
dist_noinst_SCRIPTS = xc-run_testsuite xc-reference.pl
TESTS = xc-run_testsuite

xc_get_data_SOURCES = xc-get_data.c
xc_get_data_LDADD = -L../src/ -lxc -lm
xc_get_data_CPPFLAGS = -I../src/

xc_consistency_SOURCES = xc-consistency.c
xc_consistency_LDADD = -L../src/ -lxc -lm
xc_consistency_CPPFLAGS = -I../src/

dist_noinst_DATA = *.data

CLEANFILES = *~ *.bak
DISTCLEANFILES = Makefile.in
