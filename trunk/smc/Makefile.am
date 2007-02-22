INCLUDES = -I$(top_srcdir)/detection -I$(top_srcdir)/estimation \
	-I$(top_srcdir)/randomutil -I$(top_srcdir)/src -I$(top_srcdir)/statistics -I$(top_srcdir)/util
METASOURCES = AUTO
lib_LIBRARIES = libsmc.a
noinst_HEADERS = ByChannelOrderResamplingAlgorithm.h Particle.h \
	ParticleFilter.h ParticleFilterWithChannelOrder.h ParticleWithChannelEstimation.h \
	ParticleWithChannelEstimationAndChannelOrder.h ParticleWithChannelEstimationAndChannelOrderAPP.h \
	ParticleWithChannelEstimationAndLinearDetection.h ParticleWithChannelEstimationAndLinearDetectionAndChannelOrderAPP.h \
	ParticleWithChannelEstimationAndLinearDetectionAndChannelOrderEstimation.h ResamplingAlgorithm.h ResamplingCriterion.h StdResamplingAlgorithm.h \
	WithChannelOrderAppParticleAddon.h WithChannelOrderEstimationParticleAddon.h \
	WithChannelOrderParticleAddon.h WithLinearDetectionParticleAddon.h smcExceptions.h
libsmc_a_SOURCES = ByChannelOrderResamplingAlgorithm.cpp Particle.cpp \
	ParticleFilter.cpp ParticleFilterWithChannelOrder.cpp ParticleWithChannelEstimation.cpp \
	ParticleWithChannelEstimationAndChannelOrder.cpp ParticleWithChannelEstimationAndChannelOrderAPP.cpp \
	ParticleWithChannelEstimationAndLinearDetection.cpp ParticleWithChannelEstimationAndLinearDetectionAndChannelOrderAPP.cpp \
	ParticleWithChannelEstimationAndLinearDetectionAndChannelOrderEstimation.cpp ResamplingCriterion.cpp StdResamplingAlgorithm.cpp \
	WithChannelOrderAppParticleAddon.cpp WithChannelOrderEstimationParticleAddon.cpp \
	WithChannelOrderParticleAddon.cpp WithLinearDetectionParticleAddon.cpp
