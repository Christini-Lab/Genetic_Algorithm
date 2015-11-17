CXXFLAGS = -O3 -g -fopenmp
LFLAGS =
LIBS =

TOP_DIR := ../../$(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))/
GA_SRC_DIR = ../../src/
GA_INCLUDE_DIR = ../../include/

GA_SRCFILES := $(wildcard $(GA_SRC_DIR)*.cpp)
GA_OBJFILES := $(patsubst %.cpp,%.o,$(wildcard $(GA_SRC_DIR)*.cpp))
GA_DEPS := $(wildcard $(GA_INCLUDE_DIR)*.hpp)

SRCFILES = $(GA_SRCFILES) main.cpp

EXEC = run

all: $(EXEC)
.PHONY: all clean

# Automatically determine dependencies of source files
depend: .depend
.depend: $(SRCFILES)
	$(CXX) $(CXXFLAGS) -MM $^ >> ./.depend;
-include .depend

$(EXEC): $(GA_OBJFILES) main.cpp
	$(CXX) $(CXXFLAGS) $(LFLAGS) -o $@ $^ $(LIBS)

clean:
	-rm -f *.o $(GA_OBJFILES) $(EXEC) .depend
