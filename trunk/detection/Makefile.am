INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/util -I/usr/include/eigen2/
METASOURCES = AUTO
lib_LIBRARIES = libdetection.a
noinst_HEADERS = DecorrelatorDetector.h LinearDetector.h MMSEDetector.h \
	RMMSEDetector.h
libdetection_a_SOURCES = DecorrelatorDetector.cpp LinearDetector.cpp \
	MMSEDetector.cpp RMMSEDetector.cpp
SUBDIRS = trellis
