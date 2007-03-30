INCLUDES = -I$(top_srcdir)/randomutil -I$(top_srcdir)/src \
	-I$(top_srcdir)/statistics -I$(top_srcdir)/transmision -I$(top_srcdir)/util -I/usr/local/include
METASOURCES = AUTO
lib_LIBRARIES = libtransmision.a
libtransmision_a_SOURCES = Alphabet.cpp exceptions.h Bits.cpp MIMOChannel.cpp \
		ARchannel.cpp Noise.cpp ChannelDependentNoise.cpp Modulator.cpp Demodulator.cpp \
		StillMemoryMIMOChannel.cpp SparklingMemoryMIMOChannel.cpp SparklingMemoryARMIMOChannel.cpp \
		OneChannelOrderPerTransmitAtennaMIMOChannel.cpp ARoneChannelOrderPerTransmitAtennaMIMOChannel.cpp NullNoise.cpp

noinst_HEADERS = Bits.h Alphabet.h MIMOChannel.h ARchannel.h Noise.h \
		ChannelDependentNoise.h Demodulator.h StillMemoryMIMOChannel.h SparklingMemoryMIMOChannel.h \
		SparklingMemoryARMIMOChannel.h OneChannelOrderPerTransmitAtennaMIMOChannel.h \
		ARoneChannelOrderPerTransmitAtennaMIMOChannel.h NullNoise.h
libtransmision_a_LIBADD = $(top_builddir)/util/libutil.a \
	$(top_builddir)/statistics/libstatistics.a $(top_builddir)/smc/resampling/libresampling.a /usr/lib/liblapackpp.la
