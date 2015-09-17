# Copyright (c) 2003-2008 Rudi Cilibrasi, Rulers of the RHouse
# All rights reserved.     cilibrar@cilibrar.com
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
#     # Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     # Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     # Neither the name of the RHouse nor the
#       names of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE RULERS AND CONTRIBUTORS "AS IS" AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE RULERS AND CONTRIBUTORS BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

SUBDIRS=src tests
ACLOCAL_AMFLAGS = -I m4
#INSTALL_DATA=$(INSTALL)
#INSTALLdir=${datadir}/qsearch


DISTCLEAN_TARGETS = config.status config.log libqsearch-0.9.4.tar.gz configure Makefile.in Doxyfile config.log depcomp config.guess config.sub Makefile qsearch.pc qsearch-gomp.pc.in qsearch-gomp.pc config.status autom4te.cache libtool missing aclocal.m4 install-sh

distclean: clean
	$(MAKE) -C src distclean ; rm -rf $(DISTCLEAN_TARGETS)

pkgconfigdir = $(libdir)/pkgconfig
if USING_GOMP
pkgconfig_DATA = qsearch-gomp.pc
else
pkgconfig_DATA = qsearch.pc
endif

m4datadir = $(datadir)/aclocal
m4files = qsearch.m4
m4data_DATA = $(m4files)

EXTRA_DIST = \
qsearch.m4 \
$(src) \
$(m4files) \
$(examples) \
$(doc)


man_MANS = \
doc/man/maketree.1 \
doc/man/qsearch-config.1

src = \
qsearch@GOMPSTR@.pc.in


doc = \
doc/man/maketree.1 \
doc/man/qsearch-config.1 \
THANKS

examples = \
examples/README.txt

