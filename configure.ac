
cflags_save="$CFLAGS"
AC_INIT(libqsearch,1.0.8,cilibrar@cilibrar.com)
AC_CONFIG_SRCDIR(src/qsearch-tree.c)
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE
AC_CANONICAL_HOST
AC_PROG_LIBTOOL

# to prevent extra -g -O2 being added
AM_PROG_CC_C_O

AM_MAINTAINER_MODE

PKG_CHECK_MODULES([CHECK], [check >= 0.9.4], [
  EXTCPPFLAGS="$EXTCPPFLAGS $CHECK_CFLAGS"
  LIBS="$LIBS $CHECK_LIBS"
  ], [
  echo -n
]
  )

if test "x$CHECK_LIBS" = x ; then
  HAVECHECK=0
else
  HAVECHECK=1
fi


AC_PATH_PROG(GOB2, gob2, "")
if test x$GOB2 = x ; then
  HAVEGOB2=0
else
  HAVEGOB2=1
fi

AC_ARG_ENABLE(openmp,
[  --enable-openmp  Enable OpenMP to support multicore processors
],
[case "${enableval}" in
   yes) ua_openmp=true ;;
   no) ua_openmp=false ;;
    *) AC_MSG_ERROR(bad value "${enableval}" for --enable-openmp) ;;
esac], [ua_openmp=false])

AX_OPENMP([ua_foundopenmp=true],[ua_foundopenmp=false])
if test x$ua_openmp = xmaybe ; then
  ua_openmp=$ua_foundopenmp
fi

OPENMPLIBS=`pkg-config --libs gthread-2.0`
LIBS="$LIBS $OPENMPLIBS"

if test x$ua_openmp = xtrue ; then
  if test x$ua_foundopenmp = xtrue ; then
    EXTCFLAGS="$EXTCFLAGS $OPENMP_CFLAGS"
    EXTCPPFLAGS="$EXTCPPFLAGS -DREENTRANT -DQSOPENMP_ENABLED=1"
  else
    echo "Cannot find OpenMP library, exitting."
    exit 1
  fi
else
    EXTCPPFLAGS="$EXTCPPFLAGS -DQSOPENMP_ENABLED=0"
fi
have_gsl=false
if pkg-config --cflags gsl ; then
  CF=`pkg-config --cflags gsl`
  EXTCPPFLAGS="$EXTCPPFLAGS $CF"
  LIBS=`pkg-config --libs gsl`
  have_gsl=true
fi

CFLAGS="$cflags_save"

AM_PATH_COMPLEARN(0.9.9, [
  EXTCPPFLAGS="$EXTCPPFLAGS $COMPLEARN_CFLAGS"
  CFLAGS="$COMPLEARN_CFLAGS"
  LIBS="$LIBS $COMPLEARN_LIBS"
  ], [
    AC_MSG_ERROR(must have complearn >= 0.9.9 installed)
] )

if test x$have_gsl = "xfalse" ; then
AM_PATH_GSL(0.8.0, [
  EXTCPPCFLAGS="$EXTCPPFLAGS $GSL_CFLAGS"
  LIBS="$LIBS $GSL_LIBS"
  ], [
    AC_MSG_ERROR(must have gsl >= 0.8.0 installed)
] )
fi

AM_PATH_XML2(2.4.0, [
  EXTCPPFLAGS="$EXTCPPFLAGS $XML_CPPFLAGS"
  LIBS="$LIBS $XML_LIBS"
], [
    AC_MSG_ERROR(must have libxml2 >= 2.4.0 installed)
] )

AM_PATH_GLIB_2_0(2.0.0, [
  EXTCPPFLAGS="$EXTCPPFLAGS $GLIB_CFLAGS"
  LIBS="$LIBS $GLIB_LIBS"
  ], [
    AC_MSG_ERROR(must have glib >= 2.0.0 and gobject installed)
] , [ gobject gmodule gthread ] )

AC_CHECK_HEADERS([malloc.h])
AC_CHECK_HEADERS([errno.h])

PRIVATELIBS="$LIBS"

if test x$ua_openmp = xtrue ; then
  GOMPSTR="-gomp"
  GOMPFLAGS="-fopenmp"
  cp $srcdir/qsearch.pc.in qsearch$GOMPSTR.pc.in
else
  GOMPSTR=""
  GOMPFLAGS=""
fi

if (uname) | grep -i darwin 2>/dev/null 1>/dev/null ; then
  SONAMEOPTS=""
else
  SONAMEOPTS="-Wl,-soname -Wl,libqsearch${GOMPSTR}-1.0${shrext_cmds}.0"
fi
libqsearch_major=`echo $VERSION | $AWK -F. '{print $1}'`
libqsearch_minor=`echo $VERSION | $AWK -F. '{print $2}'`
libqsearch_patch=[`echo "$VERSION" | $AWK '-F[.~]' '{print $3}'`]
libqsearch_release="$libqsearch_major.$libqsearch_minor"
libqsearch_version=`expr $libqsearch_major + $libqsearch_minor`":$libqsearch_patch:$libqsearch_minor"

AC_ARG_ENABLE(optimize,
[  --enable-optimize  Enable performance optimizations
],
[case "${enableval}" in
   yes) ua_optimize=true ;;
   no)  ua_optimize=false ;;
    *) AC_MSG_ERROR(bad value "${enableval}" for --enable-optimize) ;;
esac], [ua_optimize=true])

if test x$ua_optimize = xtrue ; then
  EXTCPPFLAGS="$EXTCPPFLAGS -DNDEBUG -DG_DISABLE_ASSERT"
  EXTCFLAGS="$EXTCFLAGS -O3 -Wall"
  echo "Enabling performance optimizations."
else
  EXTCFLAGS="$EXTCFLAGS -Wall -g"
fi

AC_ARG_ENABLE(rpath-link,
[  --enable-rpath-link  Enable rpath-link to support custom library development
],
[case "${enableval}" in
   yes) ua_rpath=true ;;
   no) ua_rpath=false ;;
    *) AC_MSG_ERROR(bad value "${enableval}" for --enable-rpath-link) ;;
esac], [ua_rpath=false])
AM_CONDITIONAL([WANT_RPATH], [test x$ua_rpath = xtrue])
AM_CONDITIONAL([FOUND_CHECK], [test x$HAVECHECK = x1])
AM_CONDITIONAL([FOUND_GOB2], [test x$HAVEGOB2 = x1])
AM_CONDITIONAL([USING_GOMP], [test x$ua_openmp = xtrue])

SET_SONAME_VERSIONS()

AC_SUBST(PRIVATELIBS)
AC_SUBST(SONAMEOPTS)
AC_SUBST(EXTCFLAGS)
AC_SUBST(EXTCPPFLAGS)
AC_SUBST(GOB2)
AC_SUBST(GOMPSTR)
AC_SUBST(GOMPFLAGS)

CFLAGS="$EXTCPPFLAGS $EXTCFLAGS"
AM_CONFIG_HEADER(src/privconf.h)
AC_OUTPUT([qsearch$GOMPSTR.pc qsearch.pc src/qsearch-config.c Makefile src/Makefile tests/Makefile ])
