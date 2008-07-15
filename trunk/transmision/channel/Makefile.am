INCLUDES = -I$(top_srcdir)/estimation -I$(top_srcdir)/randomutil \
	-I$(top_srcdir)/src -I$(top_srcdir)/statistics -I$(top_srcdir)/transmision \
	-I$(top_srcdir)/transmision/channel/delayPowerProfile -I$(top_srcdir)/util -I/usr/local/include
METASOURCES = AUTO
lib_LIBRARIES = libchannel.a
libchannel_a_SOURCES = ARchannel.cpp BesselChannel.cpp EstimatedMIMOChannel.cpp \
	MIMOChannel.cpp StillMemoryMIMOChannel.cpp TimeInvariantChannel.cpp
SUBDIRS = delayPowerProfile
noinst_HEADERS = BesselChannel.h EstimatedMIMOChannel.h TimeInvariantChannel.h
