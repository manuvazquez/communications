INCLUDES = -I$(top_srcdir)/randomutil -I$(top_srcdir)/src \
	-I$(top_srcdir)/statistics -I$(top_srcdir)/transmision -I$(top_srcdir)/util -I/usr/include/eigen2/ \
	-I/usr/local/include
METASOURCES = AUTO
lib_LIBRARIES = libdelayPowerProfile.a
noinst_HEADERS = ConstantMeanDSPowerProfile.h ContinuousPowerProfile.h \
	DelayPowerProfile.h ExponentialPowerProfile.h FlatPowerProfile.h
libdelayPowerProfile_a_SOURCES = ConstantMeanDSPowerProfile.cpp \
	ContinuousPowerProfile.cpp DelayPowerProfile.cpp ExponentialPowerProfile.cpp FlatPowerProfile.cpp
