INCLUDES = -I$(top_srcdir)/randomutil -I$(top_srcdir)/src \
	-I$(top_srcdir)/transmision -I/usr/local/include
METASOURCES = AUTO
lib_LIBRARIES = libtransmision.a
libtransmision_a_SOURCES = Alfabeto.cpp excepcionesTransmision.h Bits.cpp \
	MIMOChannel.cpp

noinst_HEADERS = Bits.h Alfabeto.h MIMOChannel.h
libtransmision_a_LIBADD = $(top_builddir)/randomutil/librandomutil.a \
	/usr/local/lib/liblapackpp.la
