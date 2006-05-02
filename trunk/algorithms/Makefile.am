INCLUDES = -I$(top_srcdir)/detection -I$(top_srcdir)/estimation \
	-I$(top_srcdir)/filters -I$(top_srcdir)/randomutil -I$(top_srcdir)/smc -I$(top_srcdir)/src \
	-I$(top_srcdir)/statistics -I$(top_srcdir)/transmision -I$(top_srcdir)/util
METASOURCES = AUTO
lib_LIBRARIES = libalgorithms.a
noinst_HEADERS = Algorithm.h UnknownChannelAlgorithm.h \
		KnownChannelOrderAlgorithm.h SMCAlgorithm.h ML_SMCAlgorithm.h LinearFilterBasedSMCAlgorithm.h
libalgorithms_a_SOURCES = Algorithm.cpp UnknownChannelAlgorithm.cpp \
		KnownChannelOrderAlgorithm.cpp SMCAlgorithm.cpp ML_SMCAlgorithm.cpp LinearFilterBasedSMCAlgorithm.cpp
libalgorithms_a_LIBADD = $(top_builddir)/transmision/libtransmision.a \
	$(top_builddir)/statistics/libstatistics.a $(top_builddir)/randomutil/librandomutil.a $(top_builddir)/smc/libsmc.a
