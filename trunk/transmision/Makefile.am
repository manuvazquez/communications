INCLUDES = -I$(top_srcdir)/randomutil -I$(top_srcdir)/src \
	-I$(top_srcdir)/statistics -I$(top_srcdir)/transmision -I$(top_srcdir)/transmision/channel \
	-I$(top_srcdir)/util
METASOURCES = AUTO
lib_LIBRARIES = libtransmision.a
libtransmision_a_SOURCES = Alphabet.cpp exceptions.h Bits.cpp Noise.cpp \
	ChannelDependentNoise.cpp Modulator.cpp Demodulator.cpp NullNoise.cpp

noinst_HEADERS = Bits.h Alphabet.h Noise.h ChannelDependentNoise.h \
	Demodulator.h NullNoise.h
SUBDIRS = channel
