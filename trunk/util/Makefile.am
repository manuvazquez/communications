


lib_LIBRARIES = libutil.a
libutil_a_SOURCES = Util.cpp
noinst_HEADERS = mylapack.h
INCLUDES = -I$(top_srcdir)/randomutil -I$(top_srcdir)/src -I/usr/local/include
libutil_a_LIBADD = $(top_builddir)/randomutil/librandomutil.a \
	/usr/local/lib/liblapackpp.la
