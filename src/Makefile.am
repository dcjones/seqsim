
SUBDIRS = yaml-cpp

bin_PROGRAMS = seqsim fit_lnorm_mix

seqsim_SOURCES = seqsim.cpp \
				 common.c         common.h \
				 read_fasta.c     read_fasta.h \
				 gtf_parse.c      gtf_parse.h \
				 str_map.c        str_map.h \
                 params.cpp       params.hpp \
                 express.cpp      express.hpp \
                 perturb.cpp      perturb.hpp \
                 generate.cpp     generate.hpp \
                 transcripts.hpp  transcripts.cpp \
                 log_norm_mix.hpp log_norm_mix.cpp


seqsim_LDADD = yaml-cpp/libyaml-cpp.la -lgsl -lgslcblas -lm -lz

fit_lnorm_mix_SOURCES = fit_lnorm_mix.c gaussian_mixture.h gaussian_mixture.c

fit_lnorm_mix_LDADD = -lgsl -lgslcblas -lm

