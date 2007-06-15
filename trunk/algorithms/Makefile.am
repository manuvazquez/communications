INCLUDES = -I$(top_srcdir)/detection -I$(top_srcdir)/detection/trellis \
	-I$(top_srcdir)/estimation -I$(top_srcdir)/filters -I$(top_srcdir)/randomutil -I$(top_srcdir)/smc \
	-I$(top_srcdir)/smc/resampling -I$(top_srcdir)/src -I$(top_srcdir)/statistics \
	-I$(top_srcdir)/transmision -I$(top_srcdir)/transmision/channel -I$(top_srcdir)/util
METASOURCES = AUTO
lib_LIBRARIES = libalgorithms.a
noinst_HEADERS = Algorithm.h ChannelOrderEstimatorSMCAlgorithm.h \
	DSISoptAlgorithm.h ISIS.h KnownChannelAlgorithm.h KnownChannelOrderAlgorithm.h \
	KnownSymbolsKalmanBasedChannelEstimator.h LinearFilterBasedMKFAlgorithm.h \
	MultipleChannelEstimatorsPerParticleSMCAlgorithm.h PSPAlgorithm.h PSPBasedSMCAlgorithm.h SMCAlgorithm.h \
	TriangularizationBasedSMCAlgorithm.h UPSPBasedSMCAlgorithm.h USIS.h UnknownChannelOrderAlgorithm.h \
	ViterbiAlgorithm.h
libalgorithms_a_SOURCES = Algorithm.cpp ChannelOrderEstimatorSMCAlgorithm.cpp \
	DSISoptAlgorithm.cpp ISIS.cpp KnownChannelAlgorithm.cpp KnownChannelOrderAlgorithm.cpp \
	KnownSymbolsKalmanBasedChannelEstimator.cpp LinearFilterBasedMKFAlgorithm.cpp LinearFilterBasedSMCAlgorithm.cpp \
	MultipleChannelEstimatorsPerParticleSMCAlgorithm.cpp PSPAlgorithm.cpp PSPBasedSMCAlgorithm.cpp SMCAlgorithm.cpp \
	TriangularizationBasedSMCAlgorithm.cpp UPSPBasedSMCAlgorithm.cpp USIS.cpp UnknownChannelAlgorithm.cpp \
	UnknownChannelOrderAlgorithm.cpp ViterbiAlgorithm.cpp
SUBDIRS = USIS2SIS
