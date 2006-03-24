INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/transmision \
	-I$(top_srcdir)/randomutil
METASOURCES = AUTO
lib_LIBRARIES = libtransmision.a
libtransmision_a_SOURCES = Alfabeto.cpp excepcionesTransmision.h Bits.cpp

noinst_HEADERS = Bits.h tipos.h Alfabeto.h
libtransmision_a_LIBADD = $(top_builddir)/randomutil/librandomutil.a
