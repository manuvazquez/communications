INCLUDES = -I$(top_srcdir)/randomutil -I$(top_srcdir)/src \
	-I$(top_srcdir)/statistics -I$(top_srcdir)/util
METASOURCES = AUTO
lib_LIBRARIES = libstatistics.a
noinst_HEADERS = ARprocess.h
libstatistics_a_SOURCES = ARprocess.cpp
libstatistics_a_LIBADD = $(top_builddir)/util/libutil.a \
	$(top_builddir)/randomutil/librandomutil.a /usr/local/lib/liblapackpp.la
