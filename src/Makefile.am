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

libxc_a_SOURCES = util.c \
	lda.c lda_x.c lda_funcs.c lda_c_vwn.c lda_c_pz.c lda_c_pw.c lda_c_amgb.c \
	gga.c gga_x.c gga_lb.c gga_c_pbe.c gga_c_lyp.c gga_c_p86.c \
	mgga.c mgga_x_tpss.c mgga_c_tpss.c \
	lca.c lca_omc.c lca_lch.c \
	functionals.c xc_f.c

noinst_HEADERS = util.h
include_HEADERS = xc.h

AM_CPPFLAGS = -I../../libstring_f

CLEANFILES = *~ *.bak
DISTCLEANFILES = Makefile.in
