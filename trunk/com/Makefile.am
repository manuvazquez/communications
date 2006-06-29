INCLUDES = -I$(top_srcdir)/algorithms -I$(top_srcdir)/detection \
	-I$(top_srcdir)/estimation -I$(top_srcdir)/filters -I$(top_srcdir)/randomutil -I$(top_srcdir)/smc \
	-I$(top_srcdir)/src -I$(top_srcdir)/statistics -I$(top_srcdir)/transmision \
	-I$(top_srcdir)/util
METASOURCES = AUTO
bin_PROGRAMS = com

com_LDADD = $(top_builddir)/estimation/libestimation.a \
	$(top_builddir)/transmision/libtransmision.a $(top_builddir)/filters/libfilters.a \
	$(top_builddir)/algorithms/libalgorithms.a $(top_builddir)/smc/libsmc.a $(top_builddir)/statistics/libstatistics.a \
	$(top_builddir)/util/libutil.a $(top_builddir)/randomutil/librandomutil.a \
	$(top_builddir)/detection/libdetection.a /usr/lib/liblapackpp.la
com_SOURCES = com.cpp
