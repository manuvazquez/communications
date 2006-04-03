INCLUDES = -I$(top_srcdir)/estimation -I$(top_srcdir)/src \
	-I$(top_srcdir)/transmision
METASOURCES = AUTO
lib_LIBRARIES = libalgorithms.a
noinst_HEADERS = Algorithm.h UnknownChannelAlgorithm.h \
	KnownChannelOrderAlgorithm.h
libalgorithms_a_SOURCES = Algorithm.cpp UnknownChannelAlgorithm.cpp \
	KnownChannelOrderAlgorithm.cpp
libalgorithms_a_LIBADD = $(top_builddir)/transmision/libtransmision.a
