INCLUDES = -I$(top_srcdir)/detection -I$(top_srcdir)/estimation \
	-I$(top_srcdir)/randomutil -I$(top_srcdir)/src -I$(top_srcdir)/statistics -I$(top_srcdir)/util
METASOURCES = AUTO
lib_LIBRARIES = libsmc.a
noinst_HEADERS = Particle.h ParticleFilter.h ParticleWithChannelEstimation.h \
	ParticleWithChannelEstimationAndActiveUsers.h ParticleWithChannelEstimationAndChannelOrderAPP.h \
	ParticleWithChannelEstimationAndLinearDetection.h \
	ParticleWithChannelEstimationAndLinearDetectionAndChannelOrderEstimation.h ResamplingCriterion.h WithActiveUsersParticleAddon.h \
	WithChannelOrderAppParticleAddon.h WithChannelOrderEstimationParticleAddon.h \
	WithLinearDetectionParticleAddon.h smcExceptions.h
libsmc_a_SOURCES = Particle.cpp ParticleFilter.cpp \
	ParticleWithChannelEstimation.cpp ParticleWithChannelEstimationAndActiveUsers.cpp \
	ParticleWithChannelEstimationAndChannelOrderAPP.cpp ParticleWithChannelEstimationAndLinearDetection.cpp \
	ParticleWithChannelEstimationAndLinearDetectionAndChannelOrderEstimation.cpp ResamplingCriterion.cpp WithActiveUsersParticleAddon.cpp \
	WithChannelOrderAppParticleAddon.cpp WithChannelOrderEstimationParticleAddon.cpp \
	WithLinearDetectionParticleAddon.cpp
SUBDIRS = resampling
