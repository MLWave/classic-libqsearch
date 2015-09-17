
if FOUND_CHECK

TESTS = check_qsearch
ACLOCAL_AMFLAGS = -I m4

check_PROGRAMS = check_qsearch
check_qsearch_SOURCES = check_qsearch.c check_neighborlist.c check_tree.c check_treemaster.c
check_qsearch_CFLAGS = @CHECK_CFLAGS@ @EXTCPPFLAGS@ @EXTCFLAGS@
check_qsearch_LDADD = @CHECK_LIBS@ -lqsearch@GOMPSTR@ -L../src

endif

