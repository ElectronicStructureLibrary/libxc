noinst_LIBRARIES = libxc.a
noinst_PROGRAMS = test.x

libxc_a_SOURCES = util.c lda_x.c \
	lda.c lda_funcs.c lda_c_vwn.c lda_c_pz.c lda_c_pw.c lda_c_amgb.c \
	gga.c gga_pbe.c gga_lb.c \
	mgga.c mgga_x_tpss.c mgga_c_tpss.c \
	xc_f.c

test_x_SOURCES = $(libxc_a_SOURCES) test.c
test_x_LDADD   = -L../liboct_parser -loct_parser @GSL_LIBS@

noinst_HEADERS = util.h xc.h

AM_CPPFLAGS = -I../liboct_parser

CLEANFILES = *~ *.bak
DISTCLEANFILES = Makefile.in

config.h:
	ln -sf ../liboct/config.h .

util.c: config.h
