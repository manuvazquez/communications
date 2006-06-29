INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/util
METASOURCES = AUTO
lib_LIBRARIES = libfilters.a
noinst_HEADERS = KalmanFilter.h
libfilters_a_SOURCES = KalmanFilter.cpp
libfilters_a_LIBADD = $(top_builddir)/util/libutil.a /usr/lib/liblapackpp.la
