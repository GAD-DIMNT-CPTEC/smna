# AC_INIT takes two parameters, the package name and a version number.
AC_INIT([inctime], [1.0], [joao.gerd@inpe.br])
AC_CONFIG_MACRO_DIRS([m4])
# AC_CONFIG_SRCDIR takes one parameter, the path and filename
# for one of the source code files.
AC_CONFIG_SRCDIR([src/inctime.f90])
AM_INIT_AUTOMAKE
AC_PROG_LIBTOOL

AC_ARG_ENABLE([debug],
    AS_HELP_STRING([--enable-debug], [Enable debugging flags [default=no]]))

AC_ARG_ENABLE([opt],
    AS_HELP_STRING([--enable-opt], [Enable optiminzation flags [default=no]]))

AC_ARG_ENABLE([dbl],
    AS_HELP_STRING([--enable-dbl], [Enable double precision flags [default=no]]))

# Compiler specific thinks
AC_ARG_ENABLE([cray],
    AS_HELP_STRING([--enable-cray[=ARG]], [Enable cray environment [default=no]]))
    dnl ENABLE-CRAY $enable_cray
    if test "x$enable_cray" != "x";then
       enable_mpi=yes
       AC_PROG_FC([ftn])
       if test "x$FC" = "x";then
          AC_MSG_ERROR([No fortran compiler found !])
       fi
       case $enable_cray in
          gnu) source config/flags.gnu ;;
          cray) source config/flags.cray;;
          *)AC_MSG_ERROR([Problem : --enable-cray=$enable_cray was set, you must choose --enable-cray=gnu or --enable-cray=cray !])
       esac
    else
       AC_PROG_FC([gfortran ifort])
       if test "x$FC" = "x";then
          AC_MSG_ERROR([No fortran compiler found !])
       fi
       case $FC in
          gfortran) source config/flags.gnu;;
          ifort) source config/flags.intel;;
          *)AC_MSG_ERROR([Problem : Found $FC compiler, but this compiler doesnot have support, use gfortran or ifort !])
       esac
    fi
# Configure flags
AC_SUBST([FCFLAGS],[$DFLT_FCFLAGS])
AC_SUBST([LDFLAGS],[$LNK_FCFLAGS])
if test "x$enable_debug" = "xyes";then
   AX_APPEND_FLAG([$DBG_FCFLAGS],[FCFLAGS])
fi
if test "x$enable_opt" = "xyes";then
   AX_APPEND_FLAG([$OPT_FCFLAGS],[FCFLAGS])
fi
if test "x$enable_dbl" = "xyes";then
   AX_APPEND_FLAG([$DBL_FCFLAGS],[FCFLAGS])
fi


# Find the compiler flag to include Fortran 90 module information from 
# another directory, and store that in the FC_MODINC variable.
AC_FC_MODULE_FLAG
if test -n "$FC_MODINC"; then
  FCFLAGS="$FCFLAGS $FC_MODINC. $FC_MODINC../include"
fi

# Find the compiler flag to write Fortran 90 module information to
# another directory, and store that in the FC_MODOUT variable.
AC_FC_MODULE_OUTPUT_FLAG
FCFLAGS="$FCFLAGS $FC_MODOUT../include"

# Find the Fortran 90 module file name extension
AC_FC_MODULE_EXTENSION

# By default, this macro turns on shared libraries if they are available,
# and also enables static libraries if they don’t conflict with the shared
# libraries.
LT_INIT

# This is a reusable macro for providing --with-libfoo functionality.
# 
# REQUIRE_LIB(name,lib,testfn,description)
# name = The complete name of the library file without the extension.
# lib = The name of the library file without the 'lib' prefix and without the extension.
# testfn = One function included in the library that can be used for a test compilation.
# description = Human readable text to be displayed if the library can't be found or
# if there's a problem during the test compilation.
AC_DEFUN([REQUIRE_LIB], [ {
 AC_ARG_WITH([$1], AC_HELP_STRING([--with-$1=<path>],[Location where $4 is
installed]),[],[with_$1=default])
 AS_IF( [test "x$with_$1" != xdefault],
 [
 LDFLAGS="$LDFLAGS -L${with_$1}/lib"
 CFLAGS="$CFLAGS -I${with_$1}/include"
 ])
 AC_CHECK_LIB($2,$3,[],
 [
 AC_MSG_ERROR([$4 was not found, try specifying --with-$1])
 ])
} ] )
# The list of libraries required by the source code that are external to
# our code.
#REQUIRE_LIB(libpcap,pcap,pcap_dump_open,[Libpcap packet capture library])
# AC_OUTPUT takes a space delimited list of where the Makefiles are to be created.
# You need to pass all directories where there is source code (including the base
# directory containing all of the source code (which won't need a path pre-pended
# to the 'Makefile' keyword).
AC_OUTPUT(Makefile libmisc/Makefile src/Makefile)


