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
## $Id$

noinst_HEADERS =  \
        lda_x.c lda_xc_1d_ehwlrg.c                         \
        gga_x_g96.c gga_x_herman.c \
	gga_x_lag.c gga_x_ak13.c gga_x_am05.c \
	gga_x_b86.c gga_k_thakkar.c gga_x_b88.c            \
	gga_x_bayesian.c gga_x_c09x.c

math2c:
	for i in $(top_srcdir)/mathematica/functionals/*.m; do \
	  functional=`echo $$i | sed 's#.*/##; s#\.m##'`; \
	  $(top_srcdir)/scripts/math2c.pl $(top_srcdir) $$functional 3; \
	done

