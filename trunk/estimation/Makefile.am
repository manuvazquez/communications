INCLUDES = -I$(top_srcdir)/filters -I$(top_srcdir)/randomutil \
	-I$(top_srcdir)/src -I$(top_srcdir)/statistics -I$(top_srcdir)/transmision \
	-I$(top_srcdir)/util
METASOURCES = AUTO
lib_LIBRARIES = libestimation.a
noinst_HEADERS = APPbasedChannelOrderEstimator.h ChannelMatrixEstimator.h \
	ChannelOrderEstimator.h KalmanEstimator.h LMSEstimator.h \
	OneChannelOrderPerTransmitAtennaWrapperEstimator.h RLSEstimator.h
libestimation_a_SOURCES = APPbasedChannelOrderEstimator.cpp \
	ChannelMatrixEstimator.cpp ChannelOrderEstimator.cpp KalmanEstimator.cpp LMSEstimator.cpp \
	OneChannelOrderPerTransmitAtennaWrapperEstimator.cpp RLSEstimator.cpp
libestimation_a_LIBADD = $(top_builddir)/util/libutil.a \
	$(top_builddir)/statistics/libstatistics.a $(top_builddir)/filters/libfilters.a \
	$(top_builddir)/randomutil/librandomutil.a /usr/lib/liblapackpp.la
