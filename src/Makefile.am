
SUBDIRS = yaml-cpp

bin_PROGRAMS = seqsim

seqsim_SOURCES = seqsim.cpp \
				 common.c      common.h \
				 fastq-parse.c fastq-parse.h \
				 gtf_parse.c   gtf_parse.h \
				 str_map.c     str_map.h \
                 params.cpp    params.hpp


seqsim_LDADD = yaml-cpp/libyaml-cpp.la -lgsl -lgslcblas -lm -lz

