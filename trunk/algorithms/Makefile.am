INCLUDES = -I$(top_srcdir)/detection -I$(top_srcdir)/detection/trellis \
	-I$(top_srcdir)/estimation -I$(top_srcdir)/filters -I$(top_srcdir)/randomutil -I$(top_srcdir)/smc \
	-I$(top_srcdir)/smc/resampling -I$(top_srcdir)/src -I$(top_srcdir)/statistics \
	-I$(top_srcdir)/transmision -I$(top_srcdir)/transmision/channel -I$(top_srcdir)/util \
	-I/usr/include/eigen2/
METASOURCES = AUTO
lib_LIBRARIES = libalgorithms.a
noinst_HEADERS = Algorithm.h CDMAunknownActiveUsersSISopt.h CMEBasedAlgorithm.h \
	CMEapplyingAlgorithm.h DSISoptAlgorithm.h ISIS.h KnownChannelAlgorithm.h \
	KnownChannelOrderAlgorithm.h KnownFlatChannelAndActiveUsersOptimalAlgorithm.h \
	KnownFlatChannelOptimalAlgorithm.h KnownSymbolsCMEapplyingAlgorithm.h \
	KnownSymbolsKalmanBasedChannelEstimatorAlgorithm.h LinearFilterBasedAlgorithm.h LinearFilterBasedCMEapplyingAlgorithm.h \
	LinearFilterBasedMKFAlgorithm.h MultipleChannelEstimatorsPerParticleSMCAlgorithm.h PSPAlgorithm.h \
	PSPBasedSMCAlgorithm.h SISoptAlgorithm.h SMCAlgorithm.h TimeVaryingChannelCMEbasedAlgorithm.h \
	TriangularizationBasedSMCAlgorithm.h USIS.h UnknownActiveUsersLinearFilterBasedSMCAlgorithm.h \
	UnknownChannelOrderAlgorithm.h ViterbiAlgorithm.h
libalgorithms_a_SOURCES = Algorithm.cpp CDMAunknownActiveUsersSISopt.cpp \
	CMEBasedAlgorithm.cpp CMEapplyingAlgorithm.cpp DSISoptAlgorithm.cpp ISIS.cpp \
	KnownChannelAlgorithm.cpp KnownChannelOrderAlgorithm.cpp \
	KnownFlatChannelAndActiveUsersOptimalAlgorithm.cpp KnownFlatChannelOptimalAlgorithm.cpp \
	KnownSymbolsCMEapplyingAlgorithm.cpp KnownSymbolsKalmanBasedChannelEstimatorAlgorithm.cpp \
	LinearFilterBasedAlgorithm.cpp LinearFilterBasedCMEapplyingAlgorithm.cpp \
	LinearFilterBasedMKFAlgorithm.cpp LinearFilterBasedSMCAlgorithm.cpp MLSDmAlgorithm.cpp \
	MLSDmFeedBackAlgorithm.cpp MultipleChannelEstimatorsPerParticleSMCAlgorithm.cpp PSPAlgorithm.cpp \
	PSPBasedSMCAlgorithm.cpp SISoptAlgorithm.cpp SMCAlgorithm.cpp \
	TimeVaryingChannelCMEbasedAlgorithm.cpp TriangularizationBasedSMCAlgorithm.cpp USIS.cpp \
	UnknownActiveUsersLinearFilterBasedSMCAlgorithm.cpp UnknownChannelAlgorithm.cpp UnknownChannelOrderAlgorithm.cpp \
	ViterbiAlgorithm.cpp
SUBDIRS = USIS2SIS
