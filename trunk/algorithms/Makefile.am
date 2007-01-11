INCLUDES = -I$(top_srcdir)/detection -I$(top_srcdir)/estimation \
	-I$(top_srcdir)/filters -I$(top_srcdir)/randomutil -I$(top_srcdir)/smc -I$(top_srcdir)/src \
	-I$(top_srcdir)/statistics -I$(top_srcdir)/transmision -I$(top_srcdir)/util
METASOURCES = AUTO
lib_LIBRARIES = libalgorithms.a
noinst_HEADERS = Algorithm.h KnownChannelOrderAlgorithm.h SMCAlgorithm.h \
	ML_SMCAlgorithm.h KnownChannelAlgorithm.h ViterbiAlgorithm.h \
	KnownSymbolsKalmanBasedChannelEstimator.h UnknownChannelOrderAlgorithm.h \
	MultipleChannelEstimatorsPerParticleSMCAlgorithm.h ISIR.h UCO_SIS.h
libalgorithms_a_SOURCES = Algorithm.cpp KnownChannelOrderAlgorithm.cpp \
		SMCAlgorithm.cpp ML_SMCAlgorithm.cpp LinearFilterBasedSMCAlgorithm.cpp \
		KnownChannelAlgorithm.cpp ViterbiAlgorithm.cpp KnownSymbolsKalmanBasedChannelEstimator.cpp \
		UnknownChannelAlgorithm.cpp UnknownChannelOrderAlgorithm.cpp \
		MultipleChannelEstimatorsPerParticleSMCAlgorithm.cpp ISIR.cpp UCO_SIS.cpp CME_UCO_SIS.cpp
libalgorithms_a_LIBADD = $(top_builddir)/transmision/libtransmision.a \
	$(top_builddir)/statistics/libstatistics.a $(top_builddir)/randomutil/librandomutil.a $(top_builddir)/smc/libsmc.a
