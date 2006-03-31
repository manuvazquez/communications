bin_PROGRAMS = comunicaciones


# set the include path found by configure
INCLUDES = -I$(top_srcdir)/. -I$(top_srcdir)/filters -I$(top_srcdir)/randomutil \
	-I$(top_srcdir)/statistics -I$(top_srcdir)/transmision -I$(top_srcdir)/util $(all_includes)

# the library search path.
comunicaciones_LDFLAGS = $(all_libraries) 
comunicaciones_SOURCES = pruebas.cpp

comunicaciones_LDADD = $(top_builddir)/transmision/libtransmision.a \
	$(top_builddir)/statistics/libstatistics.a $(top_builddir)/filters/libfilters.a $(top_builddir)/util/libutil.a \
	$(top_builddir)/randomutil/librandomutil.a /usr/local/lib/liblapackpp.la
noinst_HEADERS = types.h exceptions.h
