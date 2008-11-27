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
	lda_c_vwn.c lda_c_pz.c lda_c_pw.c lda_xc_teter93.c \
	lda_c_2d_amgb.c lda_c_2d_prm08.c \
	gga.c \
	gga_x_lg93.c gga_x_pbe.c gga_x_rpbe.c gga_x_pbea.c gga_x_mpbe.c gga_x_b86.c gga_x_b86_mgc.c \
	gga_x_b88.c gga_x_g96.c gga_x_pw86.c gga_x_pw91.c gga_x_optx.c \
	gga_x_dk87.c gga_x_ft97.c gga_x_wc.c gga_x_am05.c \
	gga_perdew.c gga_c_pbe.c gga_c_lyp.c gga_c_p86.c gga_c_pw91.c gga_c_am05.c gga_c_lm.c \
	gga_xc_lb.c gga_xc_b97.c gga_xc_edf1.c gga_xc_xlyp.c gga_xc_1w.c \
	gga_x_2d_b86_mgc.c \
	hyb_gga.c hyb_gga_xc_b3pw91.c hyb_gga_xc_b3lyp.c hyb_gga_xc_b3p86.c \
	hyb_gga_xc_o3lyp.c hyb_gga_xc_pbeh.c hyb_gga_xc_x3lyp.c hyb_gga_xc_b1wc.c \
	hyb_gga_xc_b97.c \
	mgga.c \
	mgga_x_lta.c mgga_x_tpss.c mgga_xc_vsxc.c mgga_x_m06l.c mgga_x_tau_hcth.c \
	mgga_c_tpss.c \
	lca.c lca_omc.c lca_lch.c \
	mix_func.c special_functions.c util.c functionals.c 

libxc_a_FUNC_SINGLE_SOURCES = $(libxc_a_FUNC_SOURCES:.c=_s.c)

libxc_a_SOURCES = $(libxc_a_FUNC_SOURCES)

if COMPILE_FORTRAN
  libxc_a_SOURCES += xc_f.c libxc_funcs.f90 libxc.f90
endif
if COMPILE_SINGLE
  libxc_a_SOURCES += $(libxc_a_FUNC_SINGLE_SOURCES)
if COMPILE_FORTRAN
    libxc_a_SOURCES += xc_f_s.c libxc_s.f90
endif
endif

dist_noinst_SCRIPTS = get_funcs.pl

noinst_HEADERS = \
	util.h work_lda.c \
	work_gga_x.c work_gga_becke.c work_gga_x_2d.c \
	work_mgga_x.c work_mgga_c.c \
	libxc_master.F90
include_HEADERS = xc.h xc_config.h
nodist_include_HEADERS = xc_funcs.h

if COMPILE_SINGLE
  include_HEADERS += xc_s.h
endif

AM_CPPFLAGS = $(string_f_include)

SUFFIXES = _s.c
.c_s.c:
	@CPP@ -D SINGLE_PRECISION=1 @CPPFLAGS@ $(INCLUDES) $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $< > $*_s.c

CLEANFILES = *~ *.bak *.mod *.il *.d *.pc* ifc* funcs_*.c funcs.h \
	libxc.f90 libxc_funcs.f90 xc_funcs.h *_s.c *_s.f90

lda.o: xc_funcs.h xc.h xc_s.h
xc_funcs.h: get_funcs.pl
	$(srcdir)/get_funcs.pl $(srcdir) $(top_builddir)

libxc_funcs.f90: xc_funcs.h
libxc.f90: libxc_master.F90 xc_funcs.h
	@FCCPP@                       @CPPFLAGS@ $(AM_CPPFLAGS) $(srcdir)/libxc_master.F90 > $(top_builddir)/src/libxc.f90
	@if [ "@F90_ACCEPTS_LINE_NUMBERS@" = "no" ]; then \
		grep -v "^#" $(top_builddir)/src/libxc.f90 > $(top_builddir)/src/libxc.f91; \
		mv -f $(top_builddir)/src/libxc.f91 $(top_builddir)/src/libxc.f90; \
	fi
libxc_s.f90: libxc_master.F90 xc_funcs.h
	@FCCPP@ -D SINGLE_PRECISION=1 @CPPFLAGS@ $(AM_CPPFLAGS) $(srcdir)/libxc_master.F90 > $(top_builddir)/src/libxc_s.f90
	@if [ "@F90_ACCEPTS_LINE_NUMBERS@" = "no" ]; then \
		grep -v "^#" $(top_builddir)/src/libxc_s.f90 > $(top_builddir)/src/libxc_s.f91; \
		mv -f $(top_builddir)/src/libxc_s.f91 $(top_builddir)/src/libxc_s.f90; \
	fi
