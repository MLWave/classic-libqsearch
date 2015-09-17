# Configure paths for QSEARCH
# Raph Levien 98-11-18
# stolen from Manish Singh    98-9-30
# stolen back from Frank Belew
# stolen from Manish Singh
# Shamelessly stolen from Owen Taylor

dnl AM_PATH_QSEARCH([MINIMUM-VERSION, [ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]])
dnl Test for QSEARCH, and define QSEARCH_CFLAGS and QSEARCH_LIBS
dnl
AC_DEFUN([AM_PATH_QSEARCH],
[dnl 
dnl Get the cflags and libraries from the qsearch-config script
dnl
AC_ARG_WITH(qsearch-prefix,[  --with-qsearch-prefix=PFX   Prefix where QSEARCH is installed (optional)],
            qsearch_prefix="$withval", qsearch_prefix="")
AC_ARG_WITH(qsearch-exec-prefix,[  --with-qsearch-exec-prefix=PFX Exec prefix where QSEARCH is installed (optional)],
            qsearch_exec_prefix="$withval", qsearch_exec_prefix="")
AC_ARG_ENABLE(qsearchtest, [  --disable-qsearchtest       Do not try to compile and run a test QSEARCH program],
		    , enable_qsearchtest=yes)

  if test x$qsearch_exec_prefix != x ; then
     qsearch_args="$qsearch_args --exec-prefix"
     if test x${QSEARCH_CONFIG+set} != xset ; then
        QSEARCH_CONFIG=$qsearch_exec_prefix/bin/qsearch-config
     fi
  fi
  if test x$qsearch_prefix != x ; then
     qsearch_args="$qsearch_args --prefix"
     if test x${QSEARCH_CONFIG+set} != xset ; then
        QSEARCH_CONFIG=$qsearch_prefix/bin/qsearch-config
     fi
  fi

  AC_PATH_PROG(QSEARCH_CONFIG, qsearch-config, no)
  min_qsearch_version=ifelse([$1], ,0.2.5,$1)
  AC_MSG_CHECKING(for QSEARCH - version >= $min_qsearch_version)
  no_qsearch=""
  if test "$QSEARCH_CONFIG" = "no" ; then
    no_qsearch=yes
  else
    QSEARCH_CFLAGS=`$QSEARCH_CONFIG $qsearchconf_args --cflags`
    QSEARCH_LIBS=`$QSEARCH_CONFIG $qsearchconf_args --libs`

    qsearch_major_version=`$QSEARCH_CONFIG $qsearch_args --version | \
           sed 's/\([[0-9]]*\).\([[0-9]]*\).\([[0-9]]*\)/\1/'`
    qsearch_minor_version=`$QSEARCH_CONFIG $qsearch_args --version | \
           sed 's/\([[0-9]]*\).\([[0-9]]*\).\([[0-9]]*\)/\2/'`
    qsearch_micro_version=`$QSEARCH_CONFIG $qsearch_config_args --version | \
           sed 's/\([[0-9]]*\).\([[0-9]]*\).\([[0-9]]*\)/\3/'`
    if test "x$enable_qsearchtest" = "xyes" ; then
      ac_save_CFLAGS="$CFLAGS"
      ac_save_LIBS="$LIBS"
dnl      CFLAGS="$CFLAGS $QSEARCH_CFLAGS"
dnl      LIBS="$LIBS $QSEARCH_LIBS"
dnl
dnl Now check if the installed QSEARCH is sufficiently new. (Also sanity
dnl checks the results of qsearch-config to some extent
dnl
      rm -f conf.qsearchtest
      AC_TRY_RUN([
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

char*
my_strdup (char *str)
{
  char *new_str;
  
  if (str)
    {
      new_str = malloc ((strlen (str) + 1) * sizeof(char));
      strcpy (new_str, str);
    }
  else
    new_str = NULL;
  
  return new_str;
}

int main ()
{
  int major, minor, micro;
  char *tmp_version;

  system ("touch conf.qsearchtest");

  /* HP/UX 9 (%@#!) writes to sscanf strings */
  tmp_version = my_strdup("$min_qsearch_version");
  if (sscanf(tmp_version, "%d.%d.%d", &major, &minor, &micro) != 3) {
     printf("%s, bad version string\n", "$min_qsearch_version");
     exit(1);
   }

   if (($qsearch_major_version > major) ||
      (($qsearch_major_version == major) && ($qsearch_minor_version > minor)) ||
      (($qsearch_major_version == major) && ($qsearch_minor_version == minor) && ($qsearch_micro_version >= micro)))
    {
      return 0;
    }
  else
    {
      printf("\n*** 'qsearch-config --version' returned %d.%d.%d, but the minimum version\n", $qsearch_major_version, $qsearch_minor_version, $qsearch_micro_version);
      printf("*** of QSEARCH required is %d.%d.%d. If qsearch-config is correct, then it is\n", major, minor, micro);
      printf("*** best to upgrade to the required version.\n");
      printf("*** If qsearch-config was wrong, set the environment variable QSEARCH_CONFIG\n");
      printf("*** to point to the correct copy of qsearch-config, and remove the file\n");
      printf("*** config.cache before re-running configure\n");
      return 1;
    }
}

],, no_qsearch=yes,[echo $ac_n "cross compiling; assumed OK... $ac_c"])
       CFLAGS="$ac_save_CFLAGS"
       LIBS="$ac_save_LIBS"
     fi
  fi
  if test "x$no_qsearch" = x ; then
     AC_MSG_RESULT(yes)
     ifelse([$2], , :, [$2])     
  else
     AC_MSG_RESULT(no)
     if test "$QSEARCH_CONFIG" = "no" ; then
       echo "*** The qsearch-config script installed by QSEARCH could not be found"
       echo "*** If QSEARCH was installed in PREFIX, make sure PREFIX/bin is in"
       echo "*** your path, or set the QSEARCH_CONFIG environment variable to the"
       echo "*** full path to qsearch-config."
     else
       if test -f conf.qsearchtest ; then
        :
       else
          echo "*** Could not run QSEARCH test program, checking why..."
          CFLAGS="$CFLAGS $QSEARCH_CFLAGS"
          LIBS="$LIBS $QSEARCH_LIBS"
          AC_TRY_LINK([
#include <stdio.h>
#include <qsearch/qsearch.h>
],      [ return 0; ],
        [ echo "*** The test program compiled, but did not run. This usually means"
          echo "*** that the run-time linker is not finding QSEARCH or finding the wrong"
          echo "*** version of QSEARCH. If it is not finding QSEARCH, you'll need to set your"
          echo "*** LD_LIBRARY_PATH environment variable, or edit /etc/ld.so.conf to point"
          echo "*** to the installed location  Also, make sure you have run ldconfig if that"
          echo "*** is required on your system"
	  echo "***"
          echo "*** If you have an old version installed, it is best to remove it, although"
          echo "*** you may also be able to get things to work by modifying LD_LIBRARY_PATH"],
        [ echo "*** The test program failed to compile or link. See the file config.log for the"
          echo "*** exact error that occured. This usually means QSEARCH was incorrectly installed"
          echo "*** or that you have moved QSEARCH since it was installed. In the latter case, you"
          echo "*** may want to edit the qsearch-config script: $QSEARCH_CONFIG" ])
          CFLAGS="$ac_save_CFLAGS"
          LIBS="$ac_save_LIBS"
       fi
     fi
     QSEARCH_CFLAGS=""
     QSEARCH_LIBS=""
     ifelse([$3], , :, [$3])
  fi
  AC_SUBST(QSEARCH_CFLAGS)
  AC_SUBST(QSEARCH_LIBS)
  rm -f conf.qsearchtest
])
