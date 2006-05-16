INCLUDES = -I$(top_srcdir)/detection -I$(top_srcdir)/estimation \
	-I$(top_srcdir)/src
METASOURCES = AUTO
lib_LIBRARIES = libsmc.a
noinst_HEADERS = ResamplingCriterion.h smcExceptions.h StdResamplingAlgorithm.h \
	Particle.h ParticleWithChannelEstimation.h \
	ParticleWithChannelEstimationAndLinearDetection.h
libsmc_a_SOURCES = ResamplingCriterion.cpp StdResamplingAlgorithm.cpp \
	Particle.cpp ParticleWithChannelEstimation.cpp \
	ParticleWithChannelEstimationAndLinearDetection.cpp
