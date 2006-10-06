## Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
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
## $Id: configure.ac 2422 2006-09-23 20:43:15Z marques $

AC_INIT([libxc],[0.9],[octopus-devel@tddft.org],[libxc])
AC_CONFIG_SRCDIR([src/xc.h])
AM_INIT_AUTOMAKE

AC_CONFIG_HEADERS([config.h])

# Installation prefix by default
AC_PREFIX_DEFAULT([/opt/etsf])

# Checks for programs.
AC_PROG_CC
AC_PROG_RANLIB
AC_LANG_PREPROC

# Checks for header files.
AC_HEADER_STDC
AC_FUNC_ALLOCA

# Checks for size of pointer
ACX_POINTER_SIZE

dnl build/install fortran bits unless --disable-fortran is used
AC_ARG_ENABLE([fortran],
	      AS_HELP_STRING([--disable-fortran], [don't build/install fortran components]),
	      [ac_cv_build_fortran=$enableval],
	      [ac_cv_build_fortran=yes])
AM_CONDITIONAL([ENABLE_FORTRAN], [test $ac_cv_build_fortran = yes])

if test $ac_cv_build_fortran = yes; then
  dnl try to find out what is the default FORTRAN 90 compiler
  acx_save_fcflags="${FCFLAGS}"
  AC_PROG_FC([], [Fortran 90])
  if test x"$FC" = x; then
    AC_MSG_ERROR([could not find Fortran 90 compiler])
  fi

  AC_LANG_PUSH(Fortran)
  AC_FC_SRCEXT(f90)
  FCFLAGS="${acx_save_fcflags}"

  dnl Fortran default flags
  ACX_FCFLAGS

  dnl how Fortran mangles function names
  AC_FC_WRAPPERS

  dnl Need to know the size of a Fortran integer
  ACX_FC_INTEGER_SIZE
  ACX_CC_FORTRAN_INT

  AC_LANG_POP([Fortran])
fi

AC_CONFIG_FILES([Makefile src/Makefile testsuite/Makefile])

AC_OUTPUT

