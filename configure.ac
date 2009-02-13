#
# NAME:
# configure.ac
#
# PURPOSE:
# NA
#
# REQUIREMENTS:
# NA
#
# INPUT:
# NA
#
# OUTPUT:
# NA
#
# NOTES:
# This template have only been tested on Linux systems.
#
# BUGS:
# NA
#
# AUTHOR:
# Øystein Godøy, METNO/FOU, 29.09.2006 
#
# MODIFIED:
# NA
#
# CVS_ID:
# $Id: configure.ac,v 1.1 2009-02-13 23:23:13 steingod Exp $
#

AC_PREREQ(2.59)
AC_INIT([METNO SYNOPCOL], [0.1], [o.godoy@met.no])
AC_CONFIG_SRCDIR([src/synopcol.h])
##AC_CONFIG_HEADER([src/fmcol_config.h])

# ########################################################### 
# Checks for programs.
# 
AC_PROG_CC
##AC_PROG_FC([pgf90 ifort g77 gfortran])

# ########################################################## 
# Checks for libraries, both useful and required ones
#

#
# NAME:
# FM_REQLIB
#
# PURPOSE:
# Checks for libraries that are required.
#
# AUTHOR:
# Øystein Godøy, METNO/FOU, 18.09.2007 
#
AC_DEFUN([FM_REQLIB],[
saved_CPPFLAGS="$CPPFLAGS"
saved_LDFLAGS="$LDFLAGS"
saved_LIBS="$LIBS"
AC_ARG_WITH([$1],
    AC_HELP_STRING([--with-$1=DIR],
	[the location of mandatory lib$1 files and library either as DIR or INC,LIB]),
	,
	[with_$1=yes])
case $with_$1 in
    yes)
     echo "Using system implementation of lib$1"
      LIBS="-l$1 $LIBS"
     ;;
    no)
    AC_MSG_ERROR([lib$1 is required])
     ;;
    *,*)
      addincdir="`echo $with_$1 | cut -f1 -d,`"
      addlibdir="`echo $with_$1 | cut -f2 -d,`"
      CPPFLAGS="$CPPFLAGS -I$addincdir"
      LDFLAGS="$LDFLAGS -L$addlibdir"
      LIBS="-Wl,-rpath=$addlibdir -l$1 $LIBS"
#      LIBS="-l$1 $LIBS"
      ;;
    *)
      addincdir="$with_$1/include"
      addlibdir="$with_$1/lib"
      CPPFLAGS="$CPPFLAGS -I$addincdir"
      LDFLAGS="$LDFLAGS -L$addlibdir"
      LIBS="-Wl,-rpath=$addlibdir -l$1 $LIBS"
#      LIBS="-l$1 $LIBS"
      ;;
esac
if test [ $with_$1 != no]; then
    synopcol_have_feat=m4_toupper(SYNOPCOL_HAVE_LIB$1)
    AC_CHECK_LIB([$1],[$2],
        AC_DEFINE_UNQUOTED([$synopcol_have_feat],[1],[Define if $1 is available]),
    	[CPPFLAGS="$saved_CPPFLAGS";LDFLAGS="$saved_LDFLAGS";
	LIBS="$saved_LIBS";
	AC_MSG_ERROR([Did not find lib$1, this is required to continue])])
fi
])

#
# NAME:
# FM_USELIB
#
# PURPOSE:
# Checks for libraries that could be dropped.
#
# AUTHOR:
# Øystein Godøy, METNO/FOU, 18.09.2007 
#
AC_DEFUN([FM_USELIB],[
saved_CPPFLAGS="$CPPFLAGS"
saved_LDFLAGS="$LDFLAGS"
saved_LIBS="$LIBS"
AC_ARG_WITH([$1],
    AC_HELP_STRING([--with-$1=DIR],
	[the location of optional lib$1 files and library either as DIR or INC,LIB]),
	,
	[with_$1=yes])
case $with_$1 in
    yes)
     echo "Using system implementation of lib$1"
      LIBS="-l$1 $LIBS"
     ;;
    no)
    AC_MSG_WARN([Building library with lib$1 dependent functions disabled])
     ;;
    *,*)
      addincdir="`echo $with_$1 | cut -f1 -d,`"
      addlibdir="`echo $with_$1 | cut -f2 -d,`"
      CPPFLAGS="$CPPFLAGS -I$addincdir"
      LDFLAGS="$LDFLAGS -L$addlibdir"
      LIBS="-Wl,-rpath=$addlibdir -l$1 $LIBS"
#      LIBS="-l$1 $LIBS"
      ;;
    *)
      addincdir="$with_$1/include"
      addlibdir="$with_$1/lib"
      CPPFLAGS="$CPPFLAGS -I$addincdir"
      LDFLAGS="$LDFLAGS -L$addlibdir"
      LIBS="-Wl,-rpath=$addlibdir -l$1 $LIBS"
#      LIBS="-l$1 $LIBS"
      ;;
esac
if test [ $with_$1 != no]; then
    synopcol_have_feat=m4_toupper(SYNOPCOL_HAVE_LIB$1)
    AC_CHECK_LIB([$1],[$2],
        fm_have_lib$1="yes";
	AC_DEFINE_UNQUOTED([$synopcol_have_feat],[1],[Define if $1 is available]),
    	[CPPFLAGS="$saved_CPPFLAGS";LDFLAGS="$saved_LDFLAGS";
	LIBS="$saved_LIBS";
	AC_MSG_WARN([=======================================]);
	AC_MSG_WARN([Did not find lib$1]);
	AC_MSG_WARN([Disabling lib$1 dependent functions]);
	AC_MSG_WARN([=======================================])],[$3])
fi
])

FM_REQLIB([proj],[pj_init])

FM_REQLIB([hdf5],[H5Fopen])

FM_USELIB([hlhdf],[initHlHdf])

FM_USELIB([mi],[mrfelt_],[-lg2c])

if test "x$fm_have_libmi" == "xyes"; then
LIBS="$LIBS -lg2c"
FM_REQLIB([usenwp],[getfield_],[-lg2c])
#FM_USELIB([usenwp],[getfield_],[-lg2c])
fi

FM_USELIB([osihdf5],[init_osihdf])

if test "x$fm_have_libhlhdf" == "xyes"; then
FM_REQLIB([pps_common],[readCloudType])
fi

FM_REQLIB([fmutil],[fmivec])

FM_REQLIB([fmio],[fm_readheader])

FM_REQLIB([fmcol],[read123])

# ########################################################## 
# Defines some useful variables to put in config files
#
##AC_DEFINE_UNQUOTED([SYNOPCOL_PACKAGE_NAME],[$PACKAGE_NAME])
##AC_DEFINE_UNQUOTED([SYNOPCOL_PACKAGE_VERSION],[$PACKAGE_VERSION])
##AC_DEFINE_UNQUOTED([SYNOPCOL_PACKAGE_BUGREPORT],[$PACKAGE_BUGREPORT])
##AC_DEFINE_UNQUOTED([SYNOPCOL_PACKAGE_STRING],[$PACKAGE_STRING])
##AC_DEFINE_UNQUOTED([SYNOPCOL_PACKAGE_TARNAME],[$PACKAGE_TARNAME])

# ########################################################### 
# Checks for header files.
#
AC_HEADER_STDC

# ########################################################### 
# Checks for library functions.
AC_FUNC_MALLOC

AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT

# ########################################################### 
# Print Out Summary of Configuration
#
AC_MSG_RESULT([ ])
AC_MSG_RESULT([$PACKAGE_NAME $PACKAGE_VERSION])
AC_MSG_RESULT([Configuration:])
AC_MSG_RESULT([------------------------------------------------------------])
AC_MSG_RESULT([ ])
AC_MSG_RESULT([Compiler:	$CC])
AC_MSG_RESULT([ ])
AC_MSG_RESULT([Source Code:      $srcdir])
AC_MSG_RESULT([Host System:      $host])
AC_MSG_RESULT([  ])
AC_MSG_RESULT([Prefix:           $prefix])
AC_MSG_RESULT([Binaries:         $exec_prefix])
AC_MSG_RESULT([Manuals:          $mandir])
AC_MSG_RESULT([Data:             $datadir])
AC_MSG_RESULT([System Config:    $sysconfdir])
AC_MSG_RESULT([ ])
AC_MSG_RESULT([CPPFLAGS:         $CPPFLAGS])
AC_MSG_RESULT([LDFLAGS:          $LDFLAGS])
AC_MSG_RESULT([LIBS:             $LIBS])
AC_MSG_RESULT([DEFS:             $DEFS])
AC_MSG_RESULT([  ])
AC_MSG_RESULT([------------------------------------------------------------])
AC_MSG_RESULT([ ])

