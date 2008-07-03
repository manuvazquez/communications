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
	PSPvsPSPBasedSMCSystem.h PSPvsSMCSystem.h Rev2TVT2007System.h SMCSystem.h TVT2007System.h \
	TesisComplejidadReducidaARSystem.h TesisComplejidadReducidaBesselSystem.h TesisComplejidadReducidaSystem.h \
	TesisOrdenCanalDesconocidoARSystem.h TesisOrdenCanalDesconocidoBesselSystem.h \
	TesisOrdenCanalDesconocidoSystem.h WSA08System.h
libsys_a_SOURCES = BaseSystem.cpp ChannelOrderEstimationSystem.cpp \
	Elsevier2007ARChannelSystem.cpp Elsevier2007BesselChannelSystem.cpp Elsevier2007System.cpp \
	LMSmuTestSystem.cpp PSPvsPSPBasedSMCSystem.cpp PSPvsSMCSystem.cpp Rev2TVT2007System.cpp \
	SMCSystem.cpp TVT2007System.cpp TesisComplejidadReducidaARSystem.cpp \
	TesisComplejidadReducidaBesselSystem.cpp TesisComplejidadReducidaSystem.cpp \
	TesisOrdenCanalDesconocidoARSystem.cpp TesisOrdenCanalDesconocidoBesselSystem.cpp \
	TesisOrdenCanalDesconocidoSystem.cpp WSA08System.cpp
