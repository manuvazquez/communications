INCLUDES = -I$(top_srcdir)/randomutil -I$(top_srcdir)/src \
	-I$(top_srcdir)/statistics -I$(top_srcdir)/transmision -I$(top_srcdir)/util -I/usr/local/include
METASOURCES = AUTO
lib_LIBRARIES = libchannel.a
libchannel_a_SOURCES = ARchannel.cpp \
	ARoneChannelOrderPerTransmitAtennaMIMOChannel.cpp MIMOChannel.cpp OneChannelOrderPerTransmitAtennaMIMOChannel.cpp \
	SparklingMemoryARMIMOChannel.cpp SparklingMemoryMIMOChannel.cpp StillMemoryMIMOChannel.cpp
SUBDIRS = delayPowerProfile
