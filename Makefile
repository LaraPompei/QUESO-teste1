QUESO_DIR = /home/LIBRARIES/QUESO-0.50.0/
BOOST_DIR = /usr/local/lib
GSL_DIR = /usr/local/lib

INC_PATHS += -I. -I$(QUESO_DIR)/include -I$(BOOST_DIR)/lib/boost/ -I$(GSL_DIR)/include/gsl/ 

LIBS = \
	-L$(QUESO_DIR)/lib -lqueso \
	-L$(BOOST_DIR)/lib -lboost_program_options \
	-L$(GSL_DIR)/lib -lgsl 

CXX = mpic++
CXXFLAGS += -O3 -Wall -c

default: all

.SUFFIXES: .o .C

all: teste1_gsl

clean:
	rm -f *~
	rm -f *.o

teste1_gsl: main.o likelihood.o compute.o
	$(CXX) main.o likelihood.o compute.o -o teste1_gsl $(LIBS)

%.o:%.cpp
	echo $(INC_PATHS)
	$(CXX) $(INC_PATHS) $(CXXFLAGS) $<
