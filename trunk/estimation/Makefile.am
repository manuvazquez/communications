INCLUDES = -I$(top_srcdir)/filters -I$(top_srcdir)/src -I$(top_srcdir)/util
METASOURCES = AUTO
lib_LIBRARIES = libestimation.a
noinst_HEADERS = ChannelMatrixEstimator.h KalmanEstimator.h RLSEstimator.h \
	LMSEstimator.h
libestimation_a_SOURCES = ChannelMatrixEstimator.cpp KalmanEstimator.cpp \
		RLSEstimator.cpp LMSEstimator.cpp
libestimation_a_LIBADD = $(top_builddir)/util/libutil.a \
	$(top_builddir)/statistics/libstatistics.a $(top_builddir)/filters/libfilters.a \
	$(top_builddir)/randomutil/librandomutil.a /usr/lib/liblapackpp.la
