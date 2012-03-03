CXX = g++
#CXXFLAGS = -g -Wall -pedantic 
CXXFLAGS = -Wall -pedantic -fopenmp -O3

LIBS = -lm -lboost_program_options 
OBJS := sampler.o samplib.o stats.o frames.o
#SRCS := samplib.c sampler.cc ldaframes-sampler.cc

all: ldaframes-sampler

ldaframes-sampler: ldaframes-sampler.o $(OBJS)
	${CXX} ${CXXFLAGS} -o $@ ldaframes-sampler.o ${OBJS} ${LIBS}

.cc.o:
	${CXX} ${CXXFLAGS} -c $< 

clean:
	rm ldaframes-sampler.o samplib.o stats.o sampler.o frames.o ldaframes-sampler
