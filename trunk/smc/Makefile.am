INCLUDES = -I$(top_srcdir)/detection -I$(top_srcdir)/estimation \
	-I$(top_srcdir)/randomutil -I$(top_srcdir)/src -I$(top_srcdir)/statistics -I$(top_srcdir)/util
METASOURCES = AUTO
lib_LIBRARIES = libsmc.a
noinst_HEADERS = Particle.h ParticleFilter.h ParticleWithChannelEstimation.h \
	ParticleWithChannelEstimationAndChannelOrder.h ParticleWithChannelEstimationAndChannelOrderAPP.h \
	ParticleWithChannelEstimationAndLinearDetection.h ParticleWithChannelEstimationAndLinearDetectionAndChannelOrderAPP.h \
	ParticleWithChannelEstimationAndLinearDetectionAndChannelOrderEstimation.h ResamplingCriterion.h WithChannelOrderAppParticleAddon.h \
	WithChannelOrderEstimationParticleAddon.h WithChannelOrderParticleAddon.h WithLinearDetectionParticleAddon.h \
	smcExceptions.h
libsmc_a_SOURCES = Particle.cpp ParticleFilter.cpp \
	ParticleWithChannelEstimation.cpp ParticleWithChannelEstimationAndChannelOrder.cpp \
	ParticleWithChannelEstimationAndChannelOrderAPP.cpp ParticleWithChannelEstimationAndLinearDetection.cpp \
	ParticleWithChannelEstimationAndLinearDetectionAndChannelOrderAPP.cpp \
	ParticleWithChannelEstimationAndLinearDetectionAndChannelOrderEstimation.cpp ResamplingCriterion.cpp WithChannelOrderAppParticleAddon.cpp \
	WithChannelOrderEstimationParticleAddon.cpp WithChannelOrderParticleAddon.cpp WithLinearDetectionParticleAddon.cpp
SUBDIRS = resampling
