# $Id$

ACLOCAL_AMFLAGS = -I m4

lib_LIBRARIES = libxc.a

libxc_a_SOURCES = util.c lda_x.c \
	lda.c lda_funcs.c lda_c_vwn.c lda_c_pz.c lda_c_pw.c lda_c_amgb.c \
	gga.c gga_pbe.c gga_lb.c \
	mgga.c mgga_x_tpss.c mgga_c_tpss.c \
	lca.c lca_omc.c lca_lch.c \
	xc_f.c

noinst_HEADERS = util.h
include_HEADERS = xc.h

AM_CPPFLAGS = -I../libstring_f

CLEANFILES = *~ *.bak
DISTCLEANFILES = Makefile.in
