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
## $Id: acx.m4 3881 2008-03-12 23:51:07Z xavier $
##

################################################
# Check size of a fortran integer
# ----------------------------------
AC_DEFUN([ACX_FC_INTEGER_SIZE],
[AC_MSG_CHECKING([for the size of a Fortran integer])
  AC_REQUIRE([AC_PROG_FC])
  if test -z "$FC_INTEGER_SIZE"; then
  cat >intsizetest.f90 <<EOF
program integer_size
  integer    :: i
  integer(8) :: i8

  i8 = huge(i)
  select case(i8)
  case(127_8);                 i = 1
  case(32767_8);               i = 2
  case(2147483647_8);          i = 4
  case(9223372036854775807_8); i = 8
  end select

  write(*,'(i1)') i

end program integer_size
EOF
  ac_try='$FC $FCFLAGS $LDFLAGS -o intsizetest.x intsizetest.f90 1>&AC_FD_CC'
  if AC_TRY_EVAL(ac_try); then
    ac_try=""
    ac_fcintegersize=`./intsizetest.x`;
  else
    echo "configure: failed program was:" >&AC_FD_CC
    cat intsizetest.f90 >&AC_FD_CC
    AC_MSG_WARN(failed to compile f90 program to find the size of a Fortran integer)
    ac_fcintegersize=4
  fi
  rm -f intsizetest*
  AC_DEFINE_UNQUOTED(FC_INTEGER_SIZE, ${ac_fcintegersize}, [The size of a Fortran integer])
  AC_MSG_RESULT([${ac_fcintegersize} bytes])
fi
])

################################################
# Check which C type corresponds to Fortran int
# ----------------------------------
AC_DEFUN([ACX_CC_FORTRAN_INT],
[AC_MSG_CHECKING([for which C type corresponds to Fortran integer])
  AC_REQUIRE([ACX_FC_INTEGER_SIZE])
  AC_REQUIRE([AC_PROG_CC])
  if test -z "$CC_FORTRAN_INT"; then
  cat >ccfortranint.c <<EOF
#include <stdio.h>

int main(void)
{
  if(${ac_fcintegersize} == sizeof(char))
    printf("char");
  else if(${ac_fcintegersize} == sizeof(short))
    printf("short");
  else if(${ac_fcintegersize} == sizeof(int))
    printf("int");
  else if(${ac_fcintegersize} == sizeof(long))
    printf("long");
  return 1;
}
EOF
  ac_try='$CC $CFLAGS $LDFLAGS -o ccfortranint.x ccfortranint.c 1>&AC_FD_CC'
  if AC_TRY_EVAL(ac_try); then
    ac_try=""
  else
    echo "configure: failed program was:" >&AC_FD_CC
    cat ccfortranint.c >&AC_FD_CC
    rm -f ccfortranint*
    AC_MSG_ERROR(failed to compile C program to find the C type of a Fortran integer)
  fi
  ac_ccfortranint=`./ccfortranint.x`;
  rm -f ccfortranint*
  AC_DEFINE_UNQUOTED(CC_FORTRAN_INT, ${ac_ccfortranint}, [The C type of a Fortran integer])
  AC_MSG_RESULT([${ac_ccfortranint}])
fi
])


################################################
# Check whether the compiler accepts very long lines.
# ----------------------------------
AC_DEFUN([ACX_LONG_FORTRAN_LINES],
[AC_MSG_CHECKING([whether the compiler accepts very long lines])
AC_COMPILE_IFELSE( AC_LANG_PROGRAM( [], [
write(*, *) '456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678904567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789001234567890123456789045678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890' ]), 
 [acx_long_lines_ok=yes; AC_DEFINE(LONG_LINES, 1, [compiler supports long lines])], [acx_long_lines_ok=no])
AC_SUBST([LONG_LINES], [$acx_long_lines_ok])
AC_MSG_RESULT($acx_long_lines_ok)
])


################################################
# Check whether the compiler accepts preprocessor "# line-number" lines.
# ----------------------------------
AC_DEFUN([ACX_F90_ACCEPTS_LINE_NUMBERS],
[
AC_MSG_CHECKING([whether the compiler accepts "line-number" lines cast by the preprocessor])
AC_COMPILE_IFELSE(
    AC_LANG_PROGRAM( [], [# 1]),
    [acx_f90_accepts_line_numbers_ok=yes
    AC_DEFINE(F90_ACCEPTS_LINE_NUMBERS, 1, [compiler supports line-number lines])],
    [acx_f90_accepts_line_numbers_ok=no])
AC_SUBST(F90_ACCEPTS_LINE_NUMBERS, $acx_f90_accepts_line_numbers_ok)
AC_MSG_RESULT($acx_f90_accepts_line_numbers_ok)
]
)

################################################
# Check for the presence of a given function in Fortran.
# It substitutes AC_CHECK_FUNC, since the latter
# seems to fail with some autotools versions, due to a call to some broken
# version of AC_LANG_FUNC_LINK_TRY.
AC_DEFUN([ACX_FORTRAN_CHECK_FUNC],
[
AC_MSG_CHECKING([for $1])
AC_LANG_PUSH(Fortran)dnl
AC_TRY_LINK_FUNC($1, 
[
acx_fortran_check_func=yes
AC_DEFINE_UNQUOTED(AS_TR_CPP([HAVE_$1]),1, [Define if the $1 function can be called from Fortran])], 
[
acx_fortran_check_func=no
])dnl
AC_LANG_POP(Fortran)dnl
AC_MSG_RESULT($acx_fortran_check_func)
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
  AC_ARG_VAR(FCCPP, [Fortran preprocessor. Defaults to '/lib/cpp -C -ansi'])
])
