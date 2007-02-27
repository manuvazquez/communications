INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/util
METASOURCES = AUTO
lib_LIBRARIES = libdetection.a
noinst_HEADERS = LinearDetector.h RMMSEDetector.h MMSEDetector.h
libdetection_a_SOURCES = LinearDetector.cpp RMMSEDetector.cpp MMSEDetector.cpp
libdetection_a_LIBADD = /usr/lib/liblapackpp.la
SUBDIRS = trellis
