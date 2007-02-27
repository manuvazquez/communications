INCLUDES = -I$(top_srcdir)/detection -I$(top_srcdir)/detection/trellis \
	-I$(top_srcdir)/estimation -I$(top_srcdir)/filters -I$(top_srcdir)/randomutil -I$(top_srcdir)/smc \
	-I$(top_srcdir)/src -I$(top_srcdir)/statistics -I$(top_srcdir)/transmision \
	-I$(top_srcdir)/util
METASOURCES = AUTO
lib_LIBRARIES = libalgorithms.a
noinst_HEADERS = Algorithm.h DSISoptAlgorithm.h ISIS.h KnownChannelAlgorithm.h \
	KnownChannelOrderAlgorithm.h KnownSymbolsKalmanBasedChannelEstimator.h \
	MultipleChannelEstimatorsPerParticleSMCAlgorithm.h PSPAlgorithm.h SMCAlgorithm.h USIS.h UnknownChannelOrderAlgorithm.h \
	ViterbiAlgorithm.h
libalgorithms_a_SOURCES = Algorithm.cpp DSISoptAlgorithm.cpp ISIS.cpp \
	KnownChannelAlgorithm.cpp KnownChannelOrderAlgorithm.cpp \
	KnownSymbolsKalmanBasedChannelEstimator.cpp LinearFilterBasedSMCAlgorithm.cpp \
	MultipleChannelEstimatorsPerParticleSMCAlgorithm.cpp PSPAlgorithm.cpp SMCAlgorithm.cpp USIS.cpp UnknownChannelAlgorithm.cpp \
	UnknownChannelOrderAlgorithm.cpp ViterbiAlgorithm.cpp
libalgorithms_a_LIBADD = $(top_builddir)/transmision/libtransmision.a \
	$(top_builddir)/statistics/libstatistics.a $(top_builddir)/randomutil/librandomutil.a $(top_builddir)/smc/libsmc.a
SUBDIRS = USIS2SIS
