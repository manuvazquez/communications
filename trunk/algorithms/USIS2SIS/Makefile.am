INCLUDES = -I$(top_srcdir)/algorithms -I$(top_srcdir)/algorithms/USIS2SIS \
	-I$(top_srcdir)/detection -I$(top_srcdir)/estimation -I$(top_srcdir)/filters \
	-I$(top_srcdir)/randomutil -I$(top_srcdir)/smc -I$(top_srcdir)/smc/resampling -I$(top_srcdir)/src \
	-I$(top_srcdir)/statistics -I$(top_srcdir)/transmision -I$(top_srcdir)/util
METASOURCES = AUTO
lib_LIBRARIES = libUSIS2SIS.a
noinst_HEADERS = MaximumProbabilityCriterion.h TransitionCriterion.h \
	USIS2SISAlgorithm.h UniformRelatedCriterion.h
libUSIS2SIS_a_SOURCES = MaximumProbabilityCriterion.cpp USIS2SISAlgorithm.cpp \
	UniformRelatedCriterion.cpp
libUSIS2SIS_a_LIBADD = $(top_builddir)/smc/resampling/libresampling.a \
	$(top_builddir)/smc/libsmc.a
