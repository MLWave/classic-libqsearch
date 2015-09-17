bin_PROGRAMS=maketree improvetree qsearch-config

ACLOCAL_AMFLAGS = -I m4
GLOBCPPFLAGS= -I@abs_srcdir@ -I@abs_srcdir@/qsearch

if USING_GOMP
EXELDFLAGS = -lqsearch-gomp
lib_LTLIBRARIES   = libqsearch-gomp.la
else
EXELDFLAGS = -lqsearch
lib_LTLIBRARIES   = libqsearch.la
endif

BASEHEADERS=qsearch/qsearch-neighborlist.h qsearch/qsearch-tree.h qsearch/qsearch-treemaster.h qsearch/qsearch.h qsearch/qsearch-fulltree.h qsearch/clconfig.h

if WANT_RPATH
RPATHFLAGS= -rpath @abs_srcdir@/.libs
else
RPATHFLAGS=
endif

if USING_GOMP
libqsearch_gomp_la_SOURCES=qsearch-neighborlist.c qsearch-tree.c qsearch-treemaster.c qsearch-maketree.c qsearch-tree-to-dot.c qsearch-fulltree.c qsconfig.c $(BASEHEADERS)
libqsearch_gomp_la_CPPFLAGS= @EXTCPPFLAGS@ @EXTCFLAGS@ $(GLOBCPPFLAGS)
libqsearch_gomp_la_LDFLAGS= -version-info @libqsearch_version@ @SONAMEOPTS@ @LIBS@
else
libqsearch_la_SOURCES=qsearch-neighborlist.c qsearch-tree.c qsearch-treemaster.c qsearch-maketree.c qsearch-tree-to-dot.c qsearch-fulltree.c qsconfig.c $(BASEHEADERS)
libqsearch_la_CPPFLAGS= @EXTCPPFLAGS@ @EXTCFLAGS@ $(GLOBCPPFLAGS)
libqsearch_la_LDFLAGS= -version-info @libqsearch_version@ @SONAMEOPTS@ @LIBS@
endif

improvetree_SOURCES = improvetree.c
improvetree_LDADD = $(EXELDFLAGS)
improvetree_LDFLAGS = $(RPATHFLAGS)
improvetree_CPPFLAGS = @EXTCFLAGS@ @EXTCPPFLAGS@ $(GLOBCPPFLAGS)

maketree_SOURCES = maketree.c
maketree_LDADD = $(EXELDFLAGS)
maketree_LDFLAGS = $(RPATHFLAGS)
maketree_CPPFLAGS = @EXTCFLAGS@ @EXTCPPFLAGS@ $(GLOBCPPFLAGS)

qsearch_config_SOURCES = qsearch-config.c
qsearch_config_LDADD = $(EXELDFLAGS)
qsearch_config_LDFLAGS = $(RPATHFLAGS)
qsearch_config_CPPFLAGS = @EXTCFLAGS@ @EXTCPPFLAGS@ $(GLOBCPPFLAGS)

lib_includedir=$(includedir)/qsearch-1.0/qsearch
lib_include_HEADERS= $(BASEHEADERS) $(GOBHS)
lib_nobase_includedir=$(includedir)/qsearch-1.0
lib_nobase_include_HEADERS=qsearch.h

GOBCS=qsearch-maketree.c
GOBS=qsearch-maketree.gob

GOBHS=qsearch/qsearch-maketree.h qsearch/qsearch-maketree-private.h

GOBSTAMPS=qsearch-maketree.stamp

MAINTAINERCLEANFILES = $(GOBCS) $(GOBHS) $(GOBSTAMPS)
EXTRA_DIST = $(GOBCS) $(GOBHS) $(GOBS) privconf.h.in qsearch.h

if FOUND_GOB2
$(GOBHS): $(GOBSTAMPS)

$(GOBCS): $(GOBSTAMPS)

BUILT_SOURCES = $(GOBCS) $(GOBSTAMPS)

GOBSOURCES = @GOBSOURCES@

.gob.stamp:
	$(GOB2) --always-private-header --no-touch $< -o @abs_srcdir@/qsearch ; mv -f @abs_srcdir@/qsearch/*.c . ; touch $@
endif
