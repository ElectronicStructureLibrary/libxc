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
## $Id: Makefile.am 2433 2006-09-28 11:06:54Z marques $

lib_LIBRARIES = libxc.a

libxc_a_FUNC_SOURCES = \
	lda.c lda_x.c \
	lda_c_wigner.c lda_c_rpa.c lda_c_hl.c lda_c_xalpha.c \
	lda_c_vwn.c lda_c_pz.c lda_c_pw.c lda_c_amgb.c \
	gga.c \
	gga_x_lg93.c gga_x_pbe.c gga_x_rpbe.c gga_x_b86.c gga_x_b86_mgc.c \
	gga_x_b88.c gga_x_g96.c gga_x_pw86.c gga_x_pw91.c gga_x_optx.c \
	gga_x_dk87.c gga_x_ft97.c gga_x_wc.c \
	gga_c_pbe.c gga_c_lyp.c gga_c_p86.c gga_c_pw91.c \
	gga_xc_lb.c gga_xc_hcth.c gga_xc_edf1.c gga_xc_xlyp.c \
	mgga.c mgga_x_tpss.c mgga_c_tpss.c \
	lca.c lca_omc.c lca_lch.c

if COMPILE_FORTRAN
libxc_a_SOURCES = xc_f.c libxc.f90 $(libxc_a_FUNC_SOURCES) util.c functionals.c
else
libxc_a_SOURCES = $(libxc_a_FUNC_SOURCES) util.c functionals.c
endif

dist_noinst_SCRIPTS = get_funcs.pl

noinst_HEADERS = util.h work_gga_x.c
include_HEADERS = xc.h

AM_CPPFLAGS = $(string_f_include)

CLEANFILES = *~ *.bak *.mod *.il *.d *.pc* ifc* funcs_*.c funcs.h libxc_funcs.f90 xc_funcs.h

functionals.c: xc.h
util.c: xc.h
test.c: xc.h
xc_f.c: xc.h
xc.h: xc_funcs.h
libxc.o: libxc_funcs.o
libxc_funcs.f90: xc_funcs.h
xc_funcs.h: get_funcs.pl $(libxc_a_FUNC_SOURCES)
	$(srcdir)/get_funcs.pl $(srcdir) $(top_builddir)
