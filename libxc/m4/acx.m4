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
## $Id: acx.m4 2402 2006-09-14 16:59:22Z acastro $
##

################################################
# Check size of a pointer
# ----------------------------------
AC_DEFUN([ACX_POINTER_SIZE],
[AC_MSG_CHECKING([for the size of a pointer])
  AC_REQUIRE([AC_PROG_CC])
  if test -z "$POINTER_SIZE"; then
  cat >pointertest.c <<EOF
#include <stdio.h>
void main()
{
  printf("%ld", sizeof(void *));
}
EOF
  ac_try='$CC $CFLAGS -o pointertest.x pointertest.c 1>&AC_FD_CC'
  if AC_TRY_EVAL(ac_try); then
    ac_try=""
  else
    echo "configure: failed program was:" >&AC_FD_CC
    cat pointertest.c >&AC_FD_CC
    rm -f pointertest*
    AC_MSG_ERROR(failed to compile c program to find the size of a pointer)
  fi
  ac_pointersize=`./pointertest.x`;
  rm -f pointertest*
  AC_DEFINE_UNQUOTED(POINTER_SIZE, ${ac_pointersize}, [The size of a C pointer])
  AC_MSG_RESULT([${ac_pointersize} bytes])
fi
])


################################################
# Check wether the compiler accepts very long lines.
# ----------------------------------
AC_DEFUN([ACX_LONG_FORTRAN_LINES],
[AC_MSG_CHECKING([wether the compiler accepts very long lines])
AC_COMPILE_IFELSE( AC_LANG_PROGRAM( [], [
write(*, *) '45678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890'
 ]), 
 [acx_long_lines_ok=yes; AC_DEFINE(LONG_LINES, 1, [compiler supports long lines])], [acx_long_lines_ok=no] )
AC_MSG_RESULT($acx_long_lines_ok)
])


################################################
# Check for the presence of a given function in Fortran.
# It substitutes AC_CHECK_FUNC, since the latter
# seems to fail with some autotools versions, due to a call to some broken
# version of AC_LANG_FUNC_LINK_TRY.
AC_DEFUN([ACX_FORTRAN_CHECK_FUNC],
[AC_LANG_PUSH(Fortran)dnl
AC_TRY_LINK_FUNC($1, [AC_DEFINE_UNQUOTED(AS_TR_CPP([HAVE_$1]),1,
  [Define if the $1 function can be called from Fortran])], [])dnl
AC_LANG_POP(Fortran)dnl
])


################################################
# AC_LANG_FUNC_LINK_TRY(Fortran)(FUNCTION)
# ----------------------------------
m4_define([AC_LANG_FUNC_LINK_TRY(Fortran)],
[AC_LANG_PROGRAM([], [call [$1]])])

################################################
# AC_LANG_PREPROC(Fortran)
# ---------------------------
m4_define([AC_LANG_PREPROC(Fortran)],[
  # this should not be hardwired
  if test -z "$FCCPP"; then FCCPP="/lib/cpp -C -ansi"; fi
  AC_SUBST(FCCPP)
])

################################################
# Get default FFLAGS
# ----------------------------------
# this function can certainly be improved on
AC_DEFUN([ACX_FCFLAGS],
[
AC_REQUIRE([AC_CANONICAL_HOST])

if test -z "${FCFLAGS}"; then
  case "${host}" in
  i?86*linux*)
    case "${FC}" in
    pgf90*)
      FCFLAGS="-O2 -fast -Munroll -Mnoframe -Mdalign"
      ;;
    abf90*)
      FCFLAGS="-O -YEXT_NAMES=LCS -YEXT_SFX=_"
      ;;
    ifc|ifort*)
      FCFLAGS="-u -zero -fpp1 -nbs -pc80 -pad -align -unroll -O3 -ip"
      a=`echo $host | sed "s/^i//" | sed "s/86*//"`
      if test $a > 5 ; then
         FCFLAGS="$FCFLAGS -tpp7 -xW"
      fi
    ;;
    *)
      FCFLAGS="-O"
    esac
    ;;
  x86_64*)
    dnl NAG => FCFLAGS="-colour -kind=byte -mismatch_all -abi=64 -ieee=full -O4 -Ounroll=4"
    dnl ABSOFT => FCFLAGS="-O3 -mcmodel=medium -m64 -cpu:host -YEXT_NAMES=LCS -YEXT_SFX=_ -YDEALLOC=MINE"
    dnl PGI => FCFLAGS="-fastsse -mcmodel=medium -O4 -Mdalign -Mlarge_arrays -Mscalarsse -Munroll=c:4,n:4 -Mvect=assoc,sse,cachesize:262144 -Minfo"
    ;;
  alphaev*)
    FCFLAGS="-align dcommons -fast -tune host -arch host -noautomatic"
    ;;
  powerpc-ibm*)
    FCFLAGS="-bmaxdata:0x80000000 -qmaxmem=-1 -qsuffix=f=f90 -Q -O5 -qstrict -qtune=auto -qarch=auto -qhot -qipa"
    ;;
  mips-sgi-irix*)
    FCFLAGS="-O3 -INLINE -n32 -LANG:recursive=on"
    ;;
  *)
    FCFLAGS="-O"
  esac
fi
AC_MSG_NOTICE([Using FCFLAGS="$FCFLAGS"])
])
