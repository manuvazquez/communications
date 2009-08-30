INCLUDES = -I$(top_srcdir)/algorithms -I$(top_srcdir)/algorithms/USIS2SIS \
	-I$(top_srcdir)/detection -I$(top_srcdir)/detection/trellis -I$(top_srcdir)/estimation \
	-I$(top_srcdir)/filters -I$(top_srcdir)/randomutil -I$(top_srcdir)/smc \
	-I$(top_srcdir)/smc/resampling -I$(top_srcdir)/src -I$(top_srcdir)/statistics -I$(top_srcdir)/systems \
	-I$(top_srcdir)/transmision -I$(top_srcdir)/transmision/channel \
	-I$(top_srcdir)/transmision/channel/delayPowerProfile -I$(top_srcdir)/util -I/usr/include/eigen2/
METASOURCES = AUTO
bin_PROGRAMS = com

com_LDADD = $(top_builddir)/systems/libsys.a \
	$(top_builddir)/transmision/libtransmision.a $(top_builddir)/statistics/libstatistics.a \
	$(top_builddir)/util/libutil.a $(top_builddir)/transmision/channel/libchannel.a \
	$(top_builddir)/transmision/channel/delayPowerProfile/libdelayPowerProfile.a $(top_builddir)/detection/libdetection.a \
	$(top_builddir)/estimation/libestimation.a $(top_builddir)/filters/libfilters.a \
	$(top_builddir)/algorithms/libalgorithms.a $(top_builddir)/algorithms/USIS2SIS/libUSIS2SIS.a \
	$(top_builddir)/detection/trellis/libtrellis.a $(top_builddir)/smc/libsmc.a \
	$(top_builddir)/smc/resampling/libresampling.a $(top_builddir)/randomutil/librandomutil.a
com_SOURCES = com.cpp
