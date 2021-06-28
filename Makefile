QUESO_DIR = $home/LIBRARIES/QUESO-0.50.0/
BOOST_DIR = $usr/local/lib/
GSL_DIR = $usr/local/
GRV_DIR = $
TRILINOS_DIR = $

INC_PATHS = \
	-I. \
	-I$(QUESODIR)/include \
	-I$(BOOST_DIR)/lib/boost/ \
	-I$(GSL_DIR)/include/gsl/ \

LIBS = \
	-L$(QUESO_DIR)/lib -lqueso \
	-L$(BOOST_DIR)/lib -lboost_program_options \
	-L$(GSL_DIR)/lib -lgsl 

CXX = mpic++
CXXFLAGS += -O3 -Wall -c

default: all

.SUFFIXES: .o .C

all: ex_gsl

clean:
	rm -f *~
	rm -f *.o
%	rm -f example_gsl

ex_gsl: example_main.o example_likelihood.o example_qoi.o example_compute.o
	$(CXX) example_main.o example_likelihood.o example_qoi.o \
		example_compute.o -o example_gsl $(LIBS)
%.o: %.C
	$(CXX) $(INC_PATHS) $(CXXFLAGS) $<
