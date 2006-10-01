INCLUDES = -I$(top_srcdir)/detection -I$(top_srcdir)/estimation \
	-I$(top_srcdir)/randomutil -I$(top_srcdir)/src -I$(top_srcdir)/statistics -I$(top_srcdir)/util
METASOURCES = AUTO
lib_LIBRARIES = libsmc.a
noinst_HEADERS = ResamplingCriterion.h smcExceptions.h StdResamplingAlgorithm.h \
		Particle.h ParticleWithChannelEstimation.h \
		ParticleWithChannelEstimationAndLinearDetection.h ParticleFilter.h ParticleWithChannelEstimationAndChannelOrder.h \
		ResamplingAlgorithm.h ParticleFilterWithChannelOrder.h ByChannelOrderResamplingAlgorithm.h \
		ParticleWithChannelEstimationAndChannelOrderAPP.h WithLinearDetectionParticleAddon.h WithChannelOrderAppParticleAddon.h \
		WithChannelOrderParticleAddon.h ParticleWithChannelEstimationAndLinearDetectionAndChannelOrderAPP.h
libsmc_a_SOURCES = ResamplingCriterion.cpp StdResamplingAlgorithm.cpp \
		Particle.cpp ParticleWithChannelEstimation.cpp \
		ParticleWithChannelEstimationAndLinearDetection.cpp ParticleFilter.cpp ParticleWithChannelEstimationAndChannelOrder.cpp \
		ParticleFilterWithChannelOrder.cpp ByChannelOrderResamplingAlgorithm.cpp \
		ParticleWithChannelEstimationAndChannelOrderAPP.cpp WithLinearDetectionParticleAddon.cpp \
		WithChannelOrderAppParticleAddon.cpp WithChannelOrderParticleAddon.cpp \
	ParticleWithChannelEstimationAndLinearDetectionAndChannelOrderAPP.cpp
