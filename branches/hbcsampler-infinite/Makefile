CXX = g++
CXXFLAGS = -g -Wall -pedantic
#CXXFLAGS = -Wall -pedantic -fopenmp -O4

LIBS = -lm -lboost_program_options

OBJS := sampler.o samplib.o stats.o antoniak.o frames.o distributions.o sampler_init.o sampler_data.o
#SRCS := samplib.c sampler.cc ldaframes-sampler.cc

all: ldaframes-sampler

ldaframes-sampler: ldaframes-sampler.o $(OBJS)
	${CXX} ${CXXFLAGS} -o $@ ldaframes-sampler.o ${OBJS} ${LIBS}

.cc.o:
	${CXX} ${CXXFLAGS} -c $< 

clean:
	rm ldaframes-sampler.o samplib.o stats.o antoniak.o sampler.o frames.o distributions.o ldaframes-sampler sampler_init.o sampler_data.o
