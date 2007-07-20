INCLUDES = -I$(top_srcdir)/randomutil -I$(top_srcdir)/src \
	-I$(top_srcdir)/statistics -I$(top_srcdir)/transmision -I$(top_srcdir)/transmision/channel \
	-I$(top_srcdir)/transmision/channel/delayPowerProfile -I$(top_srcdir)/util
METASOURCES = AUTO
lib_LIBRARIES = libtransmision.a
libtransmision_a_SOURCES = Alphabet.cpp Bits.cpp ChannelDependentNoise.cpp \
	Demodulator.cpp Modulator.cpp Noise.cpp NullNoise.cpp PowerProfileDependentNoise.cpp \
	TransmissionUtil.cpp exceptions.h

noinst_HEADERS = Alphabet.h Bits.h ChannelDependentNoise.h Demodulator.h \
	Noise.h NullNoise.h PowerProfileDependentNoise.h TransmissionUtil.h
SUBDIRS = channel
