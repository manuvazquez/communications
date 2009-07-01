INCLUDES = -I$(top_srcdir)/filters -I$(top_srcdir)/randomutil \
	-I$(top_srcdir)/src -I$(top_srcdir)/statistics -I$(top_srcdir)/transmision \
	-I$(top_srcdir)/transmision/channel -I$(top_srcdir)/util
METASOURCES = AUTO
lib_LIBRARIES = libestimation.a
noinst_HEADERS = APPbasedChannelOrderEstimator.h CDMAKalmanEstimator.h \
	ChannelMatrixEstimator.h ChannelOrderEstimator.h KalmanEstimator.h \
	KnownChannelChannelMatrixEstimator.h KnownSymbolsKalmanEstimator.h LMSEstimator.h NLMSEstimator.h \
	RLSEstimator.h
libestimation_a_SOURCES = APPbasedChannelOrderEstimator.cpp \
	CDMAKalmanEstimator.cpp ChannelMatrixEstimator.cpp ChannelOrderEstimator.cpp \
	KalmanEstimator.cpp KnownChannelChannelMatrixEstimator.cpp KnownSymbolsKalmanEstimator.cpp \
	LMSEstimator.cpp NLMSEstimator.cpp RLSEstimator.cpp
