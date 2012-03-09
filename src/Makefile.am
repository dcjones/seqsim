
SUBDIRS = yaml-cpp

bin_PROGRAMS = seqsim

seqsim_SOURCES = seqsim.cpp \
				 common.c         common.h \
				 fastq-parse.c    fastq-parse.h \
				 gtf_parse.c      gtf_parse.h \
				 str_map.c        str_map.h \
                 params.cpp       params.hpp \
                 express.cpp      express.hpp \
                 perturb.cpp      perturb.hpp \
                 generate.cpp     generate.hpp \
                 transcripts.hpp  transcripts.cpp \
                 log_norm_mix.hpp log_norm_mix.cpp


seqsim_LDADD = yaml-cpp/libyaml-cpp.la -lgsl -lgslcblas -lm -lz

