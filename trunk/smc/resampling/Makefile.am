INCLUDES = -I$(top_srcdir)/detection -I$(top_srcdir)/estimation \
	-I$(top_srcdir)/randomutil -I$(top_srcdir)/smc -I$(top_srcdir)/src -I$(top_srcdir)/statistics \
	-I$(top_srcdir)/util
METASOURCES = AUTO
lib_LIBRARIES = libresampling.a
noinst_HEADERS = MultinomialResamplingAlgorithm.h ResamplingAlgorithm.h \
	ResidualResamplingAlgorithm.h WithThresholdResamplingAlgorithmWrapper.h
libresampling_a_SOURCES = MultinomialResamplingAlgorithm.cpp \
		ResidualResamplingAlgorithm.cpp WithThresholdResamplingAlgorithmWrapper.cpp ResamplingAlgorithm.cpp
