INCLUDES = -I$(top_srcdir)/algorithms -I$(top_srcdir)/algorithms/USIS2SIS \
	-I$(top_srcdir)/detection -I$(top_srcdir)/detection/trellis -I$(top_srcdir)/estimation \
	-I$(top_srcdir)/filters -I$(top_srcdir)/randomutil -I$(top_srcdir)/smc \
	-I$(top_srcdir)/smc/resampling -I$(top_srcdir)/src -I$(top_srcdir)/statistics \
	-I$(top_srcdir)/transmision -I$(top_srcdir)/transmision/channel \
	-I$(top_srcdir)/transmision/channel/delayPowerProfile -I$(top_srcdir)/util
METASOURCES = AUTO
lib_LIBRARIES = libsys.a
noinst_HEADERS = BaseSystem.h ChannelOrderEstimationSystem.h \
	Elsevier2007ARChannelSystem.h Elsevier2007BesselChannelSystem.h Elsevier2007System.h \
	PSPBasedChannelOrderEstimationSystem.h PSPvsPSPBasedSMCSystem.h SMCSystem.h WSA08System.h
libsys_a_SOURCES = BaseSystem.cpp ChannelOrderEstimationSystem.cpp \
	Elsevier2007ARChannelSystem.cpp Elsevier2007BesselChannelSystem.cpp Elsevier2007System.cpp \
	PSPBasedChannelOrderEstimationSystem.cpp PSPvsPSPBasedSMCSystem.cpp SMCSystem.cpp WSA08System.cpp
